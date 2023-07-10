/* -*- c++ -*- */
/*
 * Copyright 2017 Pieter Robyns, William Thenaers.
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 * 2018: patches by wilfried.philips@wphilipe.eu for low data rate and implicit header decoding
 */

#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include <gnuradio/expj.h>
#include <liquid/liquid.h>
#include <numeric>
#include <algorithm>
#include <lora/loratap.h>
#include <lora/utilities.h>
#include "decodertmp_impl.h"
#include "tables.h"
#include <string>
#include <pwd.h>
#include <math.h>

#include "dbugr.hpp"
#define GRLORA_DEBUG

namespace gr {
    namespace lora {

        decodertmp::sptr decodertmp::make(float samp_rate, uint32_t bandwidth, uint8_t sf, bool implicit, uint8_t cr, bool crc, bool reduced_rate, bool disable_drift_correction) {
            return gnuradio::get_initial_sptr
                   (new decodertmp_impl(samp_rate, bandwidth, sf, implicit, cr, crc, reduced_rate, disable_drift_correction));
        }

        /**
         * The private constructor
         */
        decodertmp_impl::decodertmp_impl(float samp_rate, uint32_t bandwidth, uint8_t sf, bool implicit, uint8_t cr, bool crc, bool reduced_rate, bool disable_drift_correction)
            : gr::sync_block("decodertmp",
                             gr::io_signature::make(1, -1, sizeof(gr_complex)),
                             gr::io_signature::make(0, 0, 0)),
            d_pwr_queue(MAX_PWR_QUEUE_SIZE) {
            // Radio config
            d_state = gr::lora::DecoderState::DETECT;

            if (sf < 6 || sf > 13) {
                std::cerr << "[LoRa Decoder] ERROR : Spreading factor should be between 6 and 12 (inclusive)!" << std::endl
                          << "                       Other values are currently not supported." << std::endl;
                exit(1);
            }

            #ifdef GRLORA_DEBUG
                d_debug_samples.open("/tmp/grlora_debug", std::ios::out | std::ios::binary);
                d_debug.open("/tmp/grlora_debug_txt", std::ios::out);
                d_dbg.attach();
            #endif

            d_bw                 = bandwidth;
            d_implicit           = implicit;
            d_reduced_rate       = reduced_rate;
            d_phdr.cr            = cr;
            d_phdr.has_mac_crc   = crc;
            d_samples_per_second = samp_rate;
            d_payload_symbols    = 0;
            d_cfo_estimation     = 0.0f;
            d_dt                 = 1.0f / d_samples_per_second;
            d_sf                 = sf;
            d_bits_per_second    = (double)d_sf * (double)(4.0 / (4.0 + d_phdr.cr)) / (1u << d_sf) * d_bw;
            d_symbols_per_second = (double)d_bw / (1u << d_sf);
            d_period             = 1.0f / (double)d_symbols_per_second;
            d_bits_per_symbol    = (double)(d_bits_per_second    / d_symbols_per_second);
            d_samples_per_symbol = (uint32_t)(d_samples_per_second / d_symbols_per_second);
            d_delay_after_sync   = d_samples_per_symbol / 4u;
            d_number_of_bins     = (uint32_t)(1u << d_sf);
            d_number_of_bins_hdr = (uint32_t)(1u << (d_sf-2));
            d_decim_factor       = d_samples_per_symbol / d_number_of_bins;
            d_energy_threshold   = 0.0f;
            d_fine_sync = 0;
            d_enable_fine_sync = !disable_drift_correction;
            set_output_multiple(2 * d_samples_per_symbol);
            pack_number = 0;
            // add -by zkw
            fft_test_flag = 0;  
            preamble_last_save.resize(d_samples_per_symbol);
            fine_sync_samples.resize(d_decim_factor + d_samples_per_symbol);

            std::cout << "Bits (nominal) per symbol: \t"      << d_bits_per_symbol    << std::endl;
            std::cout << "Bins per symbol: \t"      << d_number_of_bins     << std::endl;
            std::cout << "Samples per symbol: \t"   << d_samples_per_symbol << std::endl;
            std::cout << "Decimation: \t\t"         << d_decim_factor       << std::endl;
            if(!d_enable_fine_sync) {
                std::cout << "Warning: clock drift correction disabled" << std::endl;
            }
            if(d_implicit) {
                std::cout << "CR: \t\t"         << (int)d_phdr.cr       << std::endl;
                std::cout << "CRC: \t\t"         << (int)d_phdr.has_mac_crc       << std::endl;
            }

            // Locally generated chirps
            build_ideal_chirps();

            // FFT decoding preparations
            d_fft.resize(d_samples_per_symbol);
            d_mult_hf.resize(d_samples_per_symbol);
            d_tmp.resize(d_number_of_bins);
            d_q  = fft_create_plan(d_samples_per_symbol, &d_mult_hf[0], &d_fft[0],     LIQUID_FFT_FORWARD, 0);
            d_qr = fft_create_plan(d_number_of_bins,     &d_tmp[0],     &d_mult_hf[0], LIQUID_FFT_BACKWARD, 0);

            // Hamming coding
            fec_scheme fs = LIQUID_FEC_HAMMING84;
            d_h48_fec = fec_create(fs, NULL);

            // Register gnuradio ports
            message_port_register_out(pmt::mp("frames"));
            message_port_register_out(pmt::mp("control"));
            
            struct  passwd *pwd;
            pwd = getpwuid(getuid());
            char a[50] = "rm -rf /home/";
            char b[10] = "/FFT/*";
            strcat(strcat(a, pwd->pw_name),b);
            system(a);
        }

        /**
         * Our virtual destructor.
         */
        decodertmp_impl::~decodertmp_impl() {
            #ifdef GRLORA_DEBUG
                if (d_debug_samples.is_open())
                    d_debug_samples.close();

                if (d_debug.is_open())
                    d_debug.close();
            #endif

            fft_destroy_plan(d_q);
            fft_destroy_plan(d_qr);
            fec_destroy(d_h48_fec);
        }

        void decodertmp_impl::build_ideal_chirps(void) {
            d_downchirp.resize(d_samples_per_symbol);
            d_upchirp.resize(d_samples_per_symbol);
            d_downchirp_ifreq.resize(d_samples_per_symbol);
            d_upchirp_ifreq.resize(d_samples_per_symbol);
            d_upchirp_ifreq_v.resize(d_samples_per_symbol*3);
            gr_complex tmp[d_samples_per_symbol*3];

            const double T       = -0.5 * d_bw * d_symbols_per_second;
            const double f0      = (d_bw / 2.0);
            const double pre_dir = 2.0 * M_PI;
            double t;
            gr_complex cmx       = gr_complex(1.0f, 1.0f);

            for (uint32_t i = 0u; i < d_samples_per_symbol; i++) {
                // Width in number of samples = samples_per_symbol
                // See https://en.wikipedia.org/wiki/Chirp#Linear
                t = d_dt * i;
                d_downchirp[i] = cmx * gr_expj(pre_dir * t * (f0 + T * t));
                d_upchirp[i]   = cmx * gr_expj(pre_dir * t * (f0 + T * t) * -1.0f);
            }

            // Store instantaneous frequency
            instantaneous_frequency(&d_downchirp[0], &d_downchirp_ifreq[0], d_samples_per_symbol);
            instantaneous_frequency(&d_upchirp[0],   &d_upchirp_ifreq[0],   d_samples_per_symbol);

            samples_to_file("/tmp/downchirp", &d_downchirp[0], d_downchirp.size(), sizeof(gr_complex));
            samples_to_file("/tmp/upchirp",   &d_upchirp[0],   d_upchirp.size(),   sizeof(gr_complex));

            // Upchirp sequence
            memcpy(tmp, &d_upchirp[0], sizeof(gr_complex) * d_samples_per_symbol);
            memcpy(tmp+d_samples_per_symbol, &d_upchirp[0], sizeof(gr_complex) * d_samples_per_symbol);
            memcpy(tmp+d_samples_per_symbol*2, &d_upchirp[0], sizeof(gr_complex) * d_samples_per_symbol);
            instantaneous_frequency(tmp, &d_upchirp_ifreq_v[0], d_samples_per_symbol*3);
        }

        void decodertmp_impl::values_to_file(const std::string path, const unsigned char *v, const uint32_t length, const uint32_t ppm) {
            std::ofstream out_file;
            out_file.open(path.c_str(), std::ios::out | std::ios::app);

            for (uint32_t i = 0u; i < length; i++) {
                std::string tmp = gr::lora::to_bin(v[i], ppm);
                out_file.write(tmp.c_str(), tmp.length());
                out_file.write(" ", 1);
            }
            out_file.write("\n", 1);

            out_file.close();
        }

        void decodertmp_impl::samples_to_file(const std::string path, const gr_complex *v, const uint32_t length, const uint32_t elem_size) {
            #ifdef GRLORA_DEBUG
                std::ofstream out_file;
                out_file.open(path.c_str(), std::ios::out | std::ios::binary);

                //for(std::vector<gr_complex>::const_iterator it = v.begin(); it != v.end(); ++it) {
                for (uint32_t i = 0u; i < length; i++) {
                    out_file.write(reinterpret_cast<const char *>(&v[i]), elem_size);
                }

                out_file.close();
            #else
                (void) path;
                (void) v;
                (void) length;
                (void) elem_size;
            #endif
        }

        void decodertmp_impl::samples_to_file_app(const std::string path, const gr_complex *v, const uint32_t length, const uint32_t elem_size) {
            #ifdef GRLORA_DEBUG
                std::ofstream out_file;
                out_file.open(path.c_str(), std::ios::out | std::ios::binary | std::ios::app);

                //for(std::vector<gr_complex>::const_iterator it = v.begin(); it != v.end(); ++it) {
                for (uint32_t i = 0u; i < length; i++) {
                    out_file.write(reinterpret_cast<const char *>(&v[i]), elem_size);
                }

                out_file.close();
            #else
                (void) path;
                (void) v;
                (void) length;
                (void) elem_size;
            #endif
        }

        void decodertmp_impl::samples_debug(const gr_complex *v, const uint32_t length) {
            #ifdef GRLORA_DEBUG
                gr_complex start_indicator(0.0f, 32.0f);
                d_debug_samples.write(reinterpret_cast<const char *>(&start_indicator), sizeof(gr_complex));

                for (uint32_t i = 1u; i < length; i++) {
                    d_debug_samples.write(reinterpret_cast<const char *>(&v[i]), sizeof(gr_complex));
                }
            #else
                (void) v;
                (void) length;
            #endif
        }

        inline void decodertmp_impl::instantaneous_frequency(const gr_complex *in_samples, float *out_ifreq, const uint32_t window) {
            if (window < 2u) {
                std::cerr << "[LoRa Decoder] WARNING : window size < 2 !" << std::endl;
                return;
            }

            /* instantaneous_phase */
            for (uint32_t i = 1u; i < window; i++) {
                const float iphase_1 = std::arg(in_samples[i - 1]);
                      float iphase_2 = std::arg(in_samples[i]);

                // Unwrapped loops from liquid_unwrap_phase
                while ( (iphase_2 - iphase_1) >  M_PI ) iphase_2 -= 2.0f*M_PI;
                while ( (iphase_2 - iphase_1) < -M_PI ) iphase_2 += 2.0f*M_PI;

                out_ifreq[i - 1] = iphase_2 - iphase_1;
            }

            // Make sure there is no strong gradient if this value is accessed by mistake
            out_ifreq[window - 1] = out_ifreq[window - 2];
        }

        inline void decodertmp_impl::instantaneous_phase(const gr_complex *in_samples, float *out_iphase, const uint32_t window) {
            out_iphase[0] = std::arg(in_samples[0]);

            for (uint32_t i = 1u; i < window; i++) {
                out_iphase[i] = std::arg(in_samples[i]);
                // = the same as atan2(imag(in_samples[i]),real(in_samples[i]));

                // Unwrapped loops from liquid_unwrap_phase
                while ( (out_iphase[i] - out_iphase[i-1]) >  M_PI ) out_iphase[i] -= 2.0f*M_PI;
                while ( (out_iphase[i] - out_iphase[i-1]) < -M_PI ) out_iphase[i] += 2.0f*M_PI;
            }
        }

        float decodertmp_impl::cross_correlate_ifreq_fast(const float *samples_ifreq, const float *ideal_chirp, const uint32_t window) {
            float result = 0;
            volk_32f_x2_dot_prod_32f(&result, samples_ifreq, ideal_chirp, window);
            return result;
        }

        float decodertmp_impl::cross_correlate_fast(const gr_complex *samples, const gr_complex *ideal_chirp, const uint32_t window) {
            gr_complex result = 0;
            volk_32fc_x2_conjugate_dot_prod_32fc(&result, samples, ideal_chirp, window);
            return abs(result);
        }

        float decodertmp_impl::cross_correlate(const gr_complex *samples_1, const gr_complex *samples_2, const uint32_t window) {
            float result = 0.0f;

            for (uint32_t i = 0u; i < window; i++) {
                result += std::real(samples_1[i] * std::conj(samples_2[i]));
            }

            result /= (float)window;

            return result;
        }

        float decodertmp_impl::cross_correlate_ifreq(const float *samples_ifreq, const std::vector<float>& ideal_chirp, const uint32_t to_idx) {
            float result = 0.0f;

            const float average   = std::accumulate(samples_ifreq  , samples_ifreq + to_idx, 0.0f) / (float)(to_idx);
            const float chirp_avg = std::accumulate(&ideal_chirp[0], &ideal_chirp[to_idx]  , 0.0f) / (float)(to_idx);
            const float sd        =   stddev(samples_ifreq   , to_idx, average)
                                    * stddev(&ideal_chirp[0] , to_idx, chirp_avg);

            for (uint32_t i = 0u; i < to_idx; i++) {
                result += (samples_ifreq[i] - average) * (ideal_chirp[i] - chirp_avg) / sd;
            }

            result /= (float)(to_idx);

            return result;
        }

        void decodertmp_impl::fine_sync(const gr_complex* in_samples, int32_t bin_idx, int32_t search_space) {
            int32_t shift_ref = (bin_idx+0.3) * d_decim_factor;
            float samples_ifreq[d_samples_per_symbol];
            float max_correlation = 0.0f;
            int32_t lag = 0;

            instantaneous_frequency(in_samples, samples_ifreq, d_samples_per_symbol);

            for(int32_t i = -search_space+1; i < search_space; i++) {
                //float c = cross_correlate_fast(in_samples, &d_upchirp_v[shift_ref+i+d_samples_per_symbol], d_samples_per_symbol);
                float c = cross_correlate_ifreq_fast(samples_ifreq, &d_upchirp_ifreq_v[shift_ref+i+d_samples_per_symbol], d_samples_per_symbol);
                if(c > max_correlation) {
                     max_correlation = c;
                     lag = i;
                 }
            }

            #ifdef GRLORA_DEBUG
                d_debug << "LAG : " << lag << std::endl;
            #endif

            d_fine_sync = -lag;

            // Soft limit impact of correction
            
            // if(lag > 0)
            //     d_fine_sync = std::min(-lag / 2, -1);
            // else if(lag < 0)
            //     d_fine_sync = std::max(-lag / 2, 1);

            // Hard limit impact of correction
            // if(abs(d_fine_sync) >= d_decim_factor / 2)
            //     d_fine_sync = 0;

            //d_fine_sync = 0;
            #ifdef GRLORA_DEBUG
                d_debug << "FINE: " << d_fine_sync << std::endl;
            #endif
        }

        float decodertmp_impl::detect_preamble_autocorr(const gr_complex *samples, const uint32_t window) {
            const gr_complex* chirp1 = samples;
            const gr_complex* chirp2 = samples + d_samples_per_symbol;
            float magsq_chirp1[window];
            float magsq_chirp2[window];
            float energy_chirp1 = 0;
            float energy_chirp2 = 0;
            float autocorr = 0;
            gr_complex dot_product;

            volk_32fc_x2_conjugate_dot_prod_32fc(&dot_product, chirp1, chirp2, window);
            volk_32fc_magnitude_squared_32f(magsq_chirp1, chirp1, window);
            volk_32fc_magnitude_squared_32f(magsq_chirp2, chirp2, window);
            volk_32f_accumulator_s32f(&energy_chirp1, magsq_chirp1, window);
            volk_32f_accumulator_s32f(&energy_chirp2, magsq_chirp2, window);

            // When using implicit mode, stop when energy is halved.
            d_energy_threshold = energy_chirp2 / 2.0f;

            // For calculating the SNR later on
            d_pwr_queue.push_back(energy_chirp1 / d_samples_per_symbol);

            // Autocorr value
            autocorr = abs(dot_product / gr_complex(sqrt(energy_chirp1 * energy_chirp2), 0));

            return autocorr;
        }

        float decodertmp_impl::determine_energy(const gr_complex *samples) {
            float magsq_chirp[d_samples_per_symbol];
            float energy_chirp = 0;
            volk_32fc_magnitude_squared_32f(magsq_chirp, samples, d_samples_per_symbol);
            volk_32f_accumulator_s32f(&energy_chirp, magsq_chirp, d_samples_per_symbol);

            return energy_chirp;
        }

        void decodertmp_impl::determine_snr() {
            if(d_pwr_queue.size() >= 2) {
                float pwr_noise = d_pwr_queue[0];
                float pwr_signal = d_pwr_queue[d_pwr_queue.size()-1];
                d_snr = pwr_signal / pwr_noise;
            }
        }

        float decodertmp_impl::detect_downchirp(const gr_complex *samples, const uint32_t window) {
            float samples_ifreq[window];
            instantaneous_frequency(samples, samples_ifreq, window);

            return cross_correlate_ifreq(samples_ifreq, d_downchirp_ifreq, window - 1u);
        }

        float decodertmp_impl::detect_upchirp(const gr_complex *samples, const uint32_t window, int32_t *index) {
            float samples_ifreq[window*2];
            instantaneous_frequency(samples, samples_ifreq, window*2);

            return sliding_norm_cross_correlate_upchirp(samples_ifreq, window, index);
        }

        float decodertmp_impl::sliding_norm_cross_correlate_upchirp(const float *samples_ifreq, const uint32_t window, int32_t *index) {
             float max_correlation = 0;

             // Cross correlate
             for (uint32_t i = 0; i < window; i++) {
                 const float max_corr = cross_correlate_ifreq_fast(samples_ifreq + i, &d_upchirp_ifreq[0], window - 1u);

                 if (max_corr > max_correlation) {
                     *index = i;
                     max_correlation = max_corr;
                 }
             }

             return max_correlation;
         }

        float decodertmp_impl::stddev(const float *values, const uint32_t len, const float mean) {
            float variance = 0.0f;

            for (uint32_t i = 0u; i < len; i++) {
                const float temp = values[i] - mean;
                variance += temp * temp;
            }

            variance /= (float)len;
            return std::sqrt(variance);
        }

        /**
         *  Currently unstable due to center frequency offset.
         */
        uint32_t decodertmp_impl::get_shift_fft(const gr_complex *samples) {
            float fft_mag[d_number_of_bins];
            std::vector<gr_complex> store_tmp; 
            store_tmp.resize(d_number_of_bins);

            //Multiply with ideal downchirp
            for (uint32_t i = 0u; i < d_samples_per_symbol; i++) {
                d_mult_hf[i] = samples[i] * d_downchirp[i];
            }

            // Perform FFT
            fft_execute(d_q);
            // Decimate. Note: assumes fft size is multiple of decimation factor and number of bins is even
            // This decimation should be identical to numpy's approach
            const uint32_t N = d_number_of_bins;
            memcpy(&d_tmp[0],               &d_fft[0],                                      N / 2u * sizeof(gr_complex));
            memcpy(&d_tmp[ N / 2u ],        &d_fft[d_samples_per_symbol - (N / 2u)],        N / 2u * sizeof(gr_complex));
            
            memcpy(&store_tmp[0],           &d_fft[d_samples_per_symbol-N],                 N / 2u * sizeof(gr_complex));
            memcpy(&store_tmp[ N / 2u ],    &d_fft[(N / 2u)],                               N / 2u * sizeof(gr_complex));
            
            for (uint32_t i = 0u; i < N; i++) {
                d_tmp[i] = d_tmp[i] + store_tmp[i];
            }

            // Get magnitude
            for (uint32_t i = 0u; i < d_number_of_bins; i++) {
                fft_mag[i] = std::abs(d_tmp[i]);
            }

            // Return argmax here
            // std::cout << " " << *std::max_element(fft_mag, fft_mag + d_number_of_bins) << " ";
            return (std::max_element(fft_mag, fft_mag + d_number_of_bins) - fft_mag);
        }

        uint32_t decodertmp_impl::max_frequency_gradient_idx(const gr_complex *samples) {
            float samples_ifreq[d_samples_per_symbol];
            float samples_ifreq_avg[d_number_of_bins];

            samples_to_file("/tmp/data", &samples[0], d_samples_per_symbol, sizeof(gr_complex));

            instantaneous_frequency(samples, samples_ifreq, d_samples_per_symbol);

            for(uint32_t i = 0; i < d_number_of_bins; i++) {
                volk_32f_accumulator_s32f(&samples_ifreq_avg[i], &samples_ifreq[i*d_decim_factor], d_decim_factor);
                samples_ifreq_avg[i] /= d_decim_factor;
            }

            float max_gradient = 0.1f;
            float gradient = 0.0f;
            uint32_t max_index = 0;
            for (uint32_t i = 1u; i < d_number_of_bins; i++) {
                gradient = samples_ifreq_avg[i - 1] - samples_ifreq_avg[i];
                if (gradient > max_gradient) {
                    max_gradient = gradient;
                    max_index = i+1;
                }
            }

            return (d_number_of_bins - max_index) % d_number_of_bins;
        }

        bool decodertmp_impl::demodulate(const gr_complex *samples, const bool is_first) {
            // DBGR_TIME_MEASUREMENT_TO_FILE("SFxx_method");
            bool reduced_rate = is_first || d_reduced_rate;

            // DBGR_START_TIME_MEASUREMENT(false, "only");

            // uint32_t bin_idx = max_frequency_gradient_idx(samples);
            // std::cout << std::endl;
            // std::cout << std::dec << bin_idx << std::endl;
            uint32_t bin_idx = get_shift_fft(samples);
            // bin_idx = (bin_idx + 3) % d_number_of_bins;
            // std::cout << std::dec << bin_idx << std::endl;
            // if(d_enable_fine_sync)
            // build_ideal_chirps();
            // std::cout << (cfo_diff - 0.7) * float(d_bw)/d_number_of_bins << std::endl;
            // rebuild_ideal_chirps((cfo_diff - 0.2) * float(d_bw)/d_number_of_bins);
            // rebuild_ideal_chirps((cfo_diff - 0.2) * 122.0f);
            rebuild_ideal_chirps((cfo_diff - 0.2f) * 122.0f);
            // fine_sync(samples, bin_idx, std::max(d_decim_factor / 4u, 2u));
            
                // fine_sync(samples, bin_idx, 3);
            // fine_sync_fft(samples, bin_idx);
            fine_sync_fft2(samples, bin_idx, 1);
            rebuild_ideal_chirps(d_cfo);

            // DBGR_INTERMEDIATE_TIME_MEASUREMENT();

            // Header has additional redundancy
            if (reduced_rate) {
                bin_idx = std::lround(bin_idx / 4.0f) % d_number_of_bins_hdr;
            }

            // Decode (actually gray encode) the bin to get the symbol value
            const uint32_t word = bin_idx ^ (bin_idx >> 1u);

            #ifdef GRLORA_DEBUG
                d_debug << gr::lora::to_bin(word, reduced_rate ? d_sf - 2u : d_sf) << " " << word << " (bin " << bin_idx << ")"  << std::endl;
            #endif
            d_words.push_back(word);

            // Look for 4+cr symbols and stop

            if (d_words.size() == (4u + (is_first ? 4u : d_phdr.cr))) {
                // Deinterleave
                deinterleave(reduced_rate ? d_sf - 2u : d_sf);

                return true; // Signal that a block is ready for decoding
            }

            return false; // We need more words in order to decode a block
        }

        /**
         *  Correct the interleaving by extracting each column of bits after rotating to the left.
         *  <br/>(The words were interleaved diagonally, by rotating we make them straight into columns.)
         */
        void decodertmp_impl::deinterleave(const uint32_t ppm) {
            const uint32_t bits_per_word = d_words.size();
            const uint32_t offset_start  = ppm - 1u;

            std::vector<uint8_t> words_deinterleaved(ppm, 0u);

            if (bits_per_word > 8u) {
                // Not sure if this can ever occur. It would imply coding rate high than 4/8 e.g. 4/9.
                std::cerr << "[LoRa Decoder] WARNING : Deinterleaver: More than 8 bits per word. uint8_t will not be sufficient!\nBytes need to be stored in intermediate array and then packed into words_deinterleaved!" << std::endl;
                exit(1);
            }

            for (uint32_t i = 0u; i < bits_per_word; i++) {
                const uint32_t word = gr::lora::rotl(d_words[i], i, ppm);

                for (uint32_t j = (1u << offset_start), x = offset_start; j; j >>= 1u, x--) {
                    words_deinterleaved[x] |= !!(word & j) << i;
                }
            }

            #ifdef GRLORA_DEBUG
                print_interleave_matrix(d_debug, d_words, ppm);
                print_vector_bin(d_debug, words_deinterleaved, "D", sizeof(uint8_t) * 8u);
            #endif

            // Add to demodulated data
            d_demodulated.insert(d_demodulated.end(), words_deinterleaved.begin(), words_deinterleaved.end());

            // Cleanup
            d_words.clear();
        }

        void decodertmp_impl::decode(const bool is_header) {
            static const uint8_t shuffle_pattern[] = {5, 0, 1, 2, 4, 3, 6, 7};

            // For determining shuffle pattern
            //if (!is_header)
            //    values_to_file("/tmp/before_deshuffle", &d_demodulated[0], d_demodulated.size(), 8);

            deshuffle(shuffle_pattern, is_header);

            // For determining whitening sequence
            //if (!is_header)
            //    values_to_file("/tmp/after_deshuffle", &d_words_deshuffled[0], d_words_deshuffled.size(), 8);
            dewhiten(is_header ? gr::lora::prng_header :
                (d_phdr.cr <=2) ? gr::lora::prng_payload_cr56 : gr::lora::prng_payload_cr78);

            //if (!is_header)
            //    values_to_file("/tmp/after_dewhiten", &d_words_dewhitened[0], d_words_dewhitened.size(), 8);

            hamming_decode(is_header);
        }

        void decodertmp_impl::msg_lora_frame(void) {
            uint32_t len = sizeof(loratap_header_t) + sizeof(loraphy_header_t) + d_payload_length;
            uint32_t offset = 0;
            uint8_t buffer[len];
            loratap_header_t loratap_header;

            memset(buffer, 0, sizeof(uint8_t) * len);
            memset(&loratap_header, 0, sizeof(loratap_header));

            loratap_header.rssi.snr = (uint8_t)(10.0f * log10(d_snr) + 0.5);

            offset = gr::lora::build_packet(buffer, offset, &loratap_header, sizeof(loratap_header_t));
            offset = gr::lora::build_packet(buffer, offset, &d_phdr, sizeof(loraphy_header_t));
            offset = gr::lora::build_packet(buffer, offset, &d_decoded[0], d_payload_length);
            if(offset != len) {
                std::cerr << "decoder_impl::msg_lora_frame: invalid write" << std::endl;
                exit(1);
            }

            pmt::pmt_t payload_blob = pmt::make_blob(buffer, sizeof(uint8_t)*len);
            message_port_pub(pmt::mp("frames"), payload_blob);
        }

        void decodertmp_impl::deshuffle(const uint8_t *shuffle_pattern, const bool is_header) {
            const uint32_t to_decode = is_header ? 5u : d_demodulated.size();
            const uint32_t len       = sizeof(shuffle_pattern) / sizeof(uint8_t);
            uint8_t result;

            for (uint32_t i = 0u; i < to_decode; i++) {
                result = 0u;

                for (uint32_t j = 0u; j < len; j++) {
                    result |= !!(d_demodulated[i] & (1u << shuffle_pattern[j])) << j;
                }

                d_words_deshuffled.push_back(result);
            }

            #ifdef GRLORA_DEBUG
                print_vector_bin(d_debug, d_words_deshuffled, "S", sizeof(uint8_t)*8);
            #endif

            // We're done with these words
            if (is_header){
                d_demodulated.erase(d_demodulated.begin(), d_demodulated.begin() + 5u);
                d_words_deshuffled.push_back(0);
            } else {
                d_demodulated.clear();
            }
        }

        void decodertmp_impl::dewhiten(const uint8_t *prng) {
            const uint32_t len = d_words_deshuffled.size();

            for (uint32_t i = 0u; i < len; i++) {
                uint8_t xor_b = d_words_deshuffled[i] ^ prng[i];
                d_words_dewhitened.push_back(xor_b);
            }

            #ifdef GRLORA_DEBUG
                print_vector_bin(d_debug, d_words_dewhitened, "W", sizeof(uint8_t) * 8);
            #endif

            d_words_deshuffled.clear();
        }

        void decodertmp_impl::hamming_decode(bool is_header) {
            switch(d_phdr.cr) {
                case 4: case 3: { // Hamming(8,4) or Hamming(7,4)
                    //hamming_decode_soft(is_header);
                    uint32_t n = ceil(d_words_dewhitened.size() * 4.0f / (4.0f + d_phdr.cr));
                    uint8_t decoded[n];

                    fec_decode(d_h48_fec, n, &d_words_dewhitened[0], decoded);
                    if(!is_header)
                        swap_nibbles(decoded, n);
                    d_decoded.assign(decoded, decoded+n);
                    break;
                }
                case 2: case 1: { // Hamming(6,4) or Hamming(5,4)
                    // TODO: Report parity error to the user
                    extract_data_only(is_header);
                    break;
                }
            }

            d_words_dewhitened.clear();
        }

        /**
         * Deprecated
         */
        void decodertmp_impl::hamming_decode_soft(bool is_header) {
            uint32_t len = d_words_dewhitened.size();
            for (uint32_t i = 0u; i < len; i += 2u) {
                const uint8_t d2 = (i + 1u < len) ? hamming_decode_soft_byte(d_words_dewhitened[i + 1u]) : 0u;
                const uint8_t d1 = hamming_decode_soft_byte(d_words_dewhitened[i]);

                if(is_header)
                    d_decoded.push_back((d1 << 4u) | d2);
                else
                    d_decoded.push_back((d2 << 4u) | d1);
            }
        }

        void decodertmp_impl::extract_data_only(bool is_header) {
            static const uint8_t data_indices[4] = {1, 2, 3, 5};
            uint32_t len = d_words_dewhitened.size();

            for (uint32_t i = 0u; i < len; i += 2u) {
                const uint8_t d2 = (i + 1u < len) ? select_bits(d_words_dewhitened[i + 1u], data_indices, 4u) & 0xFF : 0u;
                const uint8_t d1 = (select_bits(d_words_dewhitened[i], data_indices, 4u) & 0xFF);

                if(is_header)
                    d_decoded.push_back((d1 << 4u) | d2);
                else
                    d_decoded.push_back((d2 << 4u) | d1);
            }
        }

        /**
         *  Old method to determine CFO. Currently unused.
         */
        void decodertmp_impl::determine_cfo(const gr_complex *samples) {
            float iphase[d_samples_per_symbol];
            const float div = (float) d_samples_per_second / (2.0f * M_PI);

            // Determine instant phase
            instantaneous_phase(samples, iphase, d_samples_per_symbol);

            float sum = 0.0f;

            for (uint32_t i = 1u; i < d_samples_per_symbol; i++) {
                sum += (float)((iphase[i] - iphase[i - 1u]) * div);
            }

            d_cfo_estimation = sum / (float)(d_samples_per_symbol - 1u);
        }

        /**
         * New method to determine CFO.
         */
        float decodertmp_impl::experimental_determine_cfo(const gr_complex *samples, uint32_t window) {
            gr_complex mult[window];
            float mult_ifreq[window];

            volk_32fc_x2_multiply_32fc(mult, samples, &d_downchirp[0], window);
            instantaneous_frequency(mult, mult_ifreq, window);

            return mult_ifreq[256] / (2.0 * M_PI) * d_samples_per_second;
        }

        float decodertmp_impl::detect_amp_rate(const gr_complex *samples, const gr_complex *ideal_chirp) {
            std::vector<float> fft_mag; 
            std::vector<gr_complex> store_tmp; 
            fft_mag.resize(d_number_of_bins);
            store_tmp.resize(d_number_of_bins);

            //Multiply with ideal downchirp
            for (uint32_t i = 0u; i < d_samples_per_symbol; i++) {
                d_mult_hf[i] = samples[i] * ideal_chirp[i];
            }

            // Perform FFT
            fft_execute(d_q);

            // Decimate. Note: assumes fft size is multiple of decimation factor and number of bins is even
            // This decimation should be identical to numpy's approach
            const uint32_t N = d_number_of_bins;
            memcpy(&d_tmp[0],               &d_fft[0],                                      N / 2u * sizeof(gr_complex));
            memcpy(&d_tmp[ N / 2u ],        &d_fft[d_samples_per_symbol - (N / 2u)],        N / 2u * sizeof(gr_complex));
            
            memcpy(&store_tmp[0],           &d_fft[d_samples_per_symbol-N],                 N / 2u * sizeof(gr_complex));
            memcpy(&store_tmp[ N / 2u ],    &d_fft[(N / 2u)],                               N / 2u * sizeof(gr_complex));
            
            for (uint32_t i = 0u; i < N; i++) {
                d_tmp[i] = d_tmp[i] + store_tmp[i];
            }

            // Get magnitude
            for (uint32_t i = 0u; i < d_number_of_bins; i++) {
                fft_mag[i] = std::abs(d_tmp[i]);
            }

            float_t max_amp = *std::max_element(fft_mag.begin(), fft_mag.end());
            float_t sum_amp = std::accumulate(fft_mag.begin(), fft_mag.end(), 0.0);
            float_t average_amp = sum_amp / fft_mag.size();
            float_t amp_rate = max_amp/average_amp;
            // calculate median
            // std::sort(fft_mag.begin(),fft_mag.end());
            // float_t median_amp = fft_mag[d_number_of_bins/2];
            // std::cout << std::dec;
            // std::cout << "The max is " << max_amp << ",ave is " << average_amp << ",rat is " << amp_rate << std::endl;
            return amp_rate;
        }

        uint32_t decodertmp_impl::get_ref_bin(const gr_complex *samples, const gr_complex *ideal_chirp) {
            float fft_mag[8192];
            std::vector<gr_complex> d_tmp_tmp;
            std::vector<gr_complex> store_tmp; 
            std::vector<gr_complex> d_fft_tmp;              ///< Vector containing the FFT resuls.
            std::vector<gr_complex> d_mult_hf_tmp;          ///< Vector containing the FFT decimation.
            fftplan d_q_tmp;  
            d_tmp_tmp.resize(8192);
            store_tmp.resize(8192);
            d_fft_tmp.resize(65536);
            d_mult_hf_tmp.resize(65536);
            d_q_tmp  = fft_create_plan(65536, &d_mult_hf_tmp[0], &d_fft_tmp[0],     LIQUID_FFT_FORWARD, 0);
            //Multiply with ideal downchirp
            for (uint32_t i = 0u; i < d_samples_per_symbol; i++) {
                d_mult_hf_tmp[i] = samples[i] * ideal_chirp[i];
            }
            for (uint32_t i = d_samples_per_symbol; i < 65536; i++) {
                d_mult_hf_tmp[i] = 0;
            }

            // Perform FFT
            fft_execute(d_q_tmp);

            // Decimate. Note: assumes fft size is multiple of decimation factor and number of bins is even
            // This decimation should be identical to numpy's approach
            const uint32_t N = 8192;
            memcpy(&d_tmp_tmp[0],           &d_fft_tmp[0],                         N / 2u * sizeof(gr_complex));
            memcpy(&d_tmp_tmp[ N / 2u ],    &d_fft_tmp[65536 - (N / 2u)],          N / 2u * sizeof(gr_complex));
            
            memcpy(&store_tmp[0],           &d_fft_tmp[65536 - N],                 N / 2u * sizeof(gr_complex));
            memcpy(&store_tmp[ N / 2u ],    &d_fft_tmp[(N / 2u)],                  N / 2u * sizeof(gr_complex));
            
            for (uint32_t i = 0u; i < N; i++) {
                d_tmp_tmp[i] = d_tmp_tmp[i] + store_tmp[i];
            }

            // Get magnitude
            for (uint32_t i = 0u; i < N; i++) {
                fft_mag[i] = std::abs(d_tmp_tmp[i]);
            }

            return (std::max_element(fft_mag, fft_mag + N) - fft_mag);
        }

        void decodertmp_impl::get_cfo_winoff(const int32_t upchirp_peak, const int32_t downchirp_peak) {
            int32_t cfo_bin = 0;
            int32_t windows_offset_tmp = 0;
            if(upchirp_peak + downchirp_peak < 8192*0.5){
                cfo_bin = upchirp_peak + downchirp_peak;
                d_cfo = -(double)cfo_bin * d_bw / (2*8192);
                windows_offset_tmp = (downchirp_peak - upchirp_peak) / pow(2,(11-d_sf));
            }
            else if(upchirp_peak + downchirp_peak > 8192*1.5){
                cfo_bin = upchirp_peak + downchirp_peak - 8192*2;
                d_cfo = -(double)cfo_bin * d_bw / (2*8192);
                windows_offset_tmp = (downchirp_peak - upchirp_peak) / pow(2,(11-d_sf));
            }
            else{
                cfo_bin = upchirp_peak + downchirp_peak - 8192;
                d_cfo = -(double)cfo_bin * d_bw / (2*8192);
                windows_offset_tmp = (8192 - (upchirp_peak - downchirp_peak)) / pow(2,(11-d_sf));
            }
            // std::cout << "windows_offset_tmp is " << windows_offset_tmp << std::endl;
            if(windows_offset_tmp < 0)
                windows_offset = d_samples_per_symbol + windows_offset_tmp;
            else windows_offset = windows_offset_tmp;
        }

        void decodertmp_impl::rebuild_ideal_chirps(const double cfo) {
            d_downchirp.resize(d_samples_per_symbol);
            d_upchirp.resize(d_samples_per_symbol);
            d_downchirp_ifreq.resize(d_samples_per_symbol);
            d_upchirp_ifreq.resize(d_samples_per_symbol);
            d_upchirp_ifreq_v.resize(d_samples_per_symbol*3);
            gr_complex tmp[d_samples_per_symbol*3];

            const double T       = -0.5 * d_bw * d_symbols_per_second;
            const double f0      = (d_bw / 2.0);
            const double pre_dir = 2.0 * M_PI;
            double t;
            gr_complex cmx       = gr_complex(1.0f, 1.0f);
            const double f1 = f0 + cfo;
            const double f2 = f0 - cfo;

            for (uint32_t i = 0u; i < d_samples_per_symbol; i++) {
                // Width in number of samples = samples_per_symbol
                // See https://en.wikipedia.org/wiki/Chirp#Linear
                t = d_dt * i;
                d_downchirp[i] = cmx * gr_expj(pre_dir * t * (f1 + T * t));
                d_upchirp[i]   = cmx * gr_expj(pre_dir * t * (f2 + T * t) * -1.0f);
            }

            // Store instantaneous frequency
            instantaneous_frequency(&d_downchirp[0], &d_downchirp_ifreq[0], d_samples_per_symbol);
            instantaneous_frequency(&d_upchirp[0],   &d_upchirp_ifreq[0],   d_samples_per_symbol);

            samples_to_file("/tmp/downchirp", &d_downchirp[0], d_downchirp.size(), sizeof(gr_complex));
            samples_to_file("/tmp/upchirp",   &d_upchirp[0],   d_upchirp.size(),   sizeof(gr_complex));

            // Upchirp sequence
            memcpy(tmp, &d_upchirp[0], sizeof(gr_complex) * d_samples_per_symbol);
            memcpy(tmp+d_samples_per_symbol, &d_upchirp[0], sizeof(gr_complex) * d_samples_per_symbol);
            memcpy(tmp+d_samples_per_symbol*2, &d_upchirp[0], sizeof(gr_complex) * d_samples_per_symbol);
            instantaneous_frequency(tmp, &d_upchirp_ifreq_v[0], d_samples_per_symbol*3);
        }

        void decodertmp_impl::fine_sync_fft(const gr_complex *samples, int32_t bin_ref) {
            float fft_mag[d_number_of_bins * 8];
            // add -by zkw
            std::vector<gr_complex> d_tmp_tmp;
            std::vector<gr_complex> store_tmp; 
            std::vector<gr_complex> d_fft_tmp;              ///< Vector containing the FFT resuls.
            std::vector<gr_complex> d_mult_hf_tmp;          ///< Vector containing the FFT decimation.
            fftplan d_q_tmp;  
            d_tmp_tmp.resize(d_number_of_bins * 8);
            store_tmp.resize(d_number_of_bins * 8);
            d_fft_tmp.resize(d_samples_per_symbol * 8);
            d_mult_hf_tmp.resize(d_samples_per_symbol * 8);
            d_q_tmp  = fft_create_plan(d_samples_per_symbol * 8, &d_mult_hf_tmp[0], &d_fft_tmp[0], LIQUID_FFT_FORWARD, 0);

            //Multiply with ideal downchirp
            for (uint32_t i = 0u; i < d_samples_per_symbol; i++) {
                d_mult_hf_tmp[i] = samples[i] * d_downchirp[i];
            }
            for (uint32_t i = d_samples_per_symbol; i < d_samples_per_symbol * 8; i++) {
                d_mult_hf_tmp[i] = 0;
            }

            // Perform FFT
            fft_execute(d_q_tmp);

            // Decimate. Note: assumes fft size is multiple of decimation factor and number of bins is even
            // This decimation should be identical to numpy's approach
            const uint32_t N = d_number_of_bins * 8;
            memcpy(&d_tmp_tmp[0],           &d_fft_tmp[0],                         N / 2u * sizeof(gr_complex));
            memcpy(&d_tmp_tmp[ N / 2u ],    &d_fft_tmp[d_samples_per_symbol*8 - (N / 2u)],          N / 2u * sizeof(gr_complex));
            
            memcpy(&store_tmp[0],           &d_fft_tmp[d_samples_per_symbol*8 - N],                 N / 2u * sizeof(gr_complex));
            memcpy(&store_tmp[ N / 2u ],    &d_fft_tmp[(N / 2u)],                  N / 2u * sizeof(gr_complex));
            
            for (uint32_t i = 0u; i < N; i++) {
                d_tmp_tmp[i] = d_tmp_tmp[i] + store_tmp[i];
            }

            // Get magnitude
            for (uint32_t i = 0u; i < N; i++) {
                fft_mag[i] = std::abs(d_tmp_tmp[i]);
            }
            int32_t max_bin = (std::max_element(fft_mag, fft_mag + N) - fft_mag);
            bin_ref = bin_ref*8 - 3;
            int32_t lag = bin_ref - max_bin;
            std::cout << "lag is " << lag << std::endl;
            d_fine_sync = -lag;
            if(lag > 0)
                d_fine_sync = std::min(-lag / 2, -1);
            else if(lag < 0)
                d_fine_sync = std::max(-lag / 2, 1);
            if(abs(d_fine_sync) >= d_decim_factor / 2)
                d_fine_sync = 0;
            // std::cout << "bin_ref is " << bin_ref << " max_bin is " << max_bin << std::endl;
            // if(d_fine_sync < -6 || d_fine_sync > 6)
            //     d_fine_sync = 0;
            // if(d_fine_sync <= -3)
            //     d_fine_sync = -1;
            // else if(d_fine_sync >= 3)  d_fine_sync = 1;
            // else d_fine_sync = 0;
        }

        void decodertmp_impl::fine_sync_fft2(const gr_complex *samples, int32_t bin_ref, int32_t search_windows){
            float fft_mag[d_number_of_bins];
            std::vector<gr_complex> store_tmp; 
            std::vector<gr_complex> samples_tmp; 
            store_tmp.resize(d_number_of_bins);
            samples_tmp.resize(d_samples_per_symbol+d_decim_factor);
            float max_value = 0, max_win = -search_windows;
            float cfo_diff =  d_cfo/122 - (uint32_t)d_cfo/122;
            // std::cout << "cfo_diff is " << cfo_diff << std::endl;
            for(int32_t win = -search_windows; win <= search_windows; win++){
                //Multiply with ideal downchirp
                if(win < 0){
                    memcpy(&samples_tmp[0],  &fine_sync_samples[0], (-win) * sizeof(gr_complex));
                    memcpy(&samples_tmp[(-win)],  &samples[0], d_samples_per_symbol * sizeof(gr_complex));
                }else   memcpy(&samples_tmp[0],  &samples[win], d_samples_per_symbol * sizeof(gr_complex));
                for (uint32_t i = 0u; i < d_samples_per_symbol; i++) {
                    d_mult_hf[i] = samples_tmp[i] * d_downchirp[i];
                }
                fft_execute(d_q);

                const uint32_t N = d_number_of_bins;
                memcpy(&d_tmp[0],               &d_fft[0],                                      N / 2u * sizeof(gr_complex));
                memcpy(&d_tmp[ N / 2u ],        &d_fft[d_samples_per_symbol - (N / 2u)],        N / 2u * sizeof(gr_complex));
                
                memcpy(&store_tmp[0],           &d_fft[d_samples_per_symbol-N],                 N / 2u * sizeof(gr_complex));
                memcpy(&store_tmp[ N / 2u ],    &d_fft[(N / 2u)],                               N / 2u * sizeof(gr_complex));
                
                for (uint32_t i = 0u; i < N; i++) {
                    d_tmp[i] = d_tmp[i] + store_tmp[i];
                }
                // Get magnitude
                for (uint32_t i = 0u; i < N; i++) {
                    fft_mag[i] = std::abs(d_tmp[i]);
                }
                float value_tmp = *std::max_element(fft_mag, fft_mag + d_number_of_bins);
                if(value_tmp >= max_value){
                    max_value = value_tmp;
                    max_win = win;
                }
            }
            d_fine_sync = max_win;
        }

        int decodertmp_impl::work(int noutput_items,
                               gr_vector_const_void_star& input_items,
                               gr_vector_void_star&       output_items) {
            (void) noutput_items;
            (void) output_items;

            const gr_complex *input     = (gr_complex *) input_items[0];
            //const gr_complex *raw_input = (gr_complex *) input_items[1]; // Input bypassed by low pass filter

            d_fine_sync = 0; // Always reset fine sync

            switch (d_state) {
                case gr::lora::DecoderState::DETECT: {
                    // float correlation = detect_preamble_autocorr(input, d_samples_per_symbol);
	  	            // std::cout << "The detect cor is "<< correlation << std::endl;
                    amp_rate = detect_amp_rate(input, &d_downchirp[0]); 
                    if (amp_rate >= 5.00f) {
                        std::cout << "Detect pass! The amp_rate is " << std::dec << amp_rate << std::endl;
                        upchirp_ref = get_ref_bin(input, &d_downchirp[0]);
                        std::cout << std::dec << "upchirp_ref is " << upchirp_ref << std::endl;
                        uint32_t upchirp_ref_tmp = get_shift_fft(input);
                        std::cout << std::dec << "upchirp_ref_tmp is " << upchirp_ref_tmp << std::endl;
                        determine_snr();
                        #ifdef GRLORA_DEBUG
                            d_debug << "Ca: " << amp_rate << std::endl;
                        #endif
                        d_corr_fails = 0u;
                        d_state = gr::lora::DecoderState::SYNC;
                        // std::cout << "DETECT pass" << std::endl;
                        consume_each(d_samples_per_symbol - upchirp_ref/pow(2,10-d_sf));

                        break;
                    }
                    consume_each(d_samples_per_symbol);

                    break;
                }

                case gr::lora::DecoderState::SYNC: {
                    // upchirp_ref = get_ref_bin(input, &d_downchirp[0]);
                    // std::cout << "upchirp_ref is " << upchirp_ref << std::endl;
                    d_state = gr::lora::DecoderState::FIND_SFD;
                    break;
                }

                case gr::lora::DecoderState::FIND_SFD: {
                    float amp_rate_SFD = detect_amp_rate(input, &d_upchirp[0]); 
                    // std::cout << "SFD amp_rate is " << amp_rate_SFD << std::endl;
                    if (amp_rate_SFD > 5.0f) {
                        upchirp_ref = get_ref_bin(&preamble_last_save[0], &d_downchirp[0]);
                        std::cout << "upchirp_ref is " << upchirp_ref << std::endl;
                        if(upchirp_ref < 128*pow(2,10-d_sf)) // sync word 0x12
                            upchirp_ref = d_samples_per_symbol*pow(2,10-d_sf) + (upchirp_ref - 128*pow(2,10-d_sf));
                        else
                            upchirp_ref = upchirp_ref - 128*pow(2,10-d_sf);
                        std::cout << "upchirp_ref is " << upchirp_ref << std::endl;
                        consume_each(d_samples_per_symbol);
                        std::cout << "SYNC pass!" <<std::endl;
                        d_state = gr::lora::DecoderState::PAUSE;
                    } else {
                        memcpy(&preamble_last_save[0],  &input[0 + d_fine_sync], d_samples_per_symbol * sizeof(gr_complex));  // save last preamble -add by zkw
                        d_corr_fails++;
                        fine_sync_fft2(input, 0, 1);
                        
                        if (d_corr_fails > 14u) {
                            d_state = gr::lora::DecoderState::DETECT;
                        }
                    }

                    consume_each(d_samples_per_symbol + d_fine_sync);
                    break;
                }

                case gr::lora::DecoderState::PAUSE: {
                    downchirp_ref = get_ref_bin(input, &d_upchirp[0]);
                    std::cout << std::dec << "The upchirp_ref is " << upchirp_ref << " The downchirp_ref is " << downchirp_ref << std::endl;
                    get_cfo_winoff(upchirp_ref, downchirp_ref);
                    std::cout << std::dec << "d_cfo is " << d_cfo << " windows_off is " << windows_offset << std::endl;
                    rebuild_ideal_chirps(d_cfo);
                    if(windows_offset < d_samples_per_symbol * 0.75){
                        memcpy(&fine_sync_samples[0],  &input[d_delay_after_sync + windows_offset - d_decim_factor], d_decim_factor*sizeof(gr_complex));
                        consume_each(d_delay_after_sync + windows_offset - d_decim_factor);
                    }
                    else{
                        memcpy(&fine_sync_samples[0],  &input[d_delay_after_sync - (d_samples_per_symbol - windows_offset) - d_decim_factor], d_decim_factor*sizeof(gr_complex));
                        consume_each(d_delay_after_sync - (d_samples_per_symbol - windows_offset) - d_decim_factor);
                    }
                    // cfo_diff = d_cfo/(float(d_bw)/d_number_of_bins) - uint32_t(d_cfo/(float(d_bw)/d_number_of_bins));
                    cfo_diff = d_cfo/122.0f - uint32_t(d_cfo/122.0f);
                        
                    std::cout << "cfo_diff is " << cfo_diff << std::endl;
                    d_state = gr::lora::DecoderState::DECODE_HEADER;
                    
                    break;
                }

                case gr::lora::DecoderState::DECODE_HEADER: {
                    if (demodulate(input, true)) {
                        if (d_implicit) {
                            d_payload_symbols = 1;
                        } else {
                            decode(true);
                            gr::lora::print_vector_hex(std::cout, &d_decoded[0], d_decoded.size(), false, false);
                            memcpy(&d_phdr, &d_decoded[0], sizeof(loraphy_header_t));
                            if (d_phdr.cr > 4)
                                d_phdr.cr = 4;
                            d_decoded.clear();

                            d_payload_length = d_phdr.length + MAC_CRC_SIZE * d_phdr.has_mac_crc;
                            //d_phy_crc = SM(decoded[1], 4, 0xf0) | MS(decoded[2], 0xf0, 4);

                            // Calculate number of payload symbols needed
                            uint8_t redundancy = (d_reduced_rate ? 2 : 0);
                            const int symbols_per_block = d_phdr.cr + 4u;
                            const float bits_needed     = float(d_payload_length) * 8.0f;
                            const float symbols_needed  = bits_needed * (symbols_per_block / 4.0f) / float(d_sf - redundancy);
                            const int blocks_needed     = (int)std::ceil(symbols_needed / symbols_per_block);
                            d_payload_symbols     = blocks_needed * symbols_per_block;

                            #ifdef GRLORA_DEBUG
                                d_debug << "LEN: " << d_payload_length << " (" << d_payload_symbols << " symbols)" << std::endl;
                            #endif
                        }
                        std::cout << "header deocded!" << std::endl;
                        d_state = gr::lora::DecoderState::DECODE_PAYLOAD;
                    }
                    // std::cout << "d_sync_sync is " << std::dec << (int32_t)d_samples_per_symbol+d_fine_sync - (int32_t)d_samples_per_symbol<< std::endl;
                    memcpy(&fine_sync_samples[0],  &input[d_samples_per_symbol+d_fine_sync-d_decim_factor], d_decim_factor*sizeof(gr_complex));
                    consume_each((int32_t)d_samples_per_symbol+d_fine_sync);
                    break;
                }

                case gr::lora::DecoderState::DECODE_PAYLOAD: {
                    if (d_implicit && determine_energy(input) < d_energy_threshold) {
                        d_payload_symbols = 0;
                        //d_demodulated.erase(d_demodulated.begin(), d_demodulated.begin() + 7u); // Test for SF 8 with header
                        d_payload_length = (int32_t)(d_demodulated.size() / 2);
                    } else if (demodulate(input, false)) {
                        if(!d_implicit)
                            d_payload_symbols -= (4u + d_phdr.cr);
                    }

                    if (d_payload_symbols <= 0) {
                        decode(false);
                        
                        gr::lora::print_vector_hex(std::cout, &d_decoded[0], d_payload_length, true, true);
                        msg_lora_frame();

                        build_ideal_chirps();
                        d_fine_sync = 0;
                        d_state = gr::lora::DecoderState::DETECT;
                        d_decoded.clear();
                        d_words.clear();
                        d_words_dewhitened.clear();
                        d_words_deshuffled.clear();
                        d_demodulated.clear();
                    }
                    // std::cout << "d_sync_sync is " << std::dec << (int32_t)d_samples_per_symbol+d_fine_sync - (int32_t)d_samples_per_symbol<< std::endl;
                    memcpy(&fine_sync_samples[0],  &input[d_samples_per_symbol+d_fine_sync-d_decim_factor], d_decim_factor*sizeof(gr_complex));
                    consume_each((int32_t)d_samples_per_symbol+d_fine_sync);
                    break;
                }

                case gr::lora::DecoderState::STOP: {
                    consume_each(d_samples_per_symbol);
                    break;
                }

                default: {
                    std::cerr << "[LoRa Decoder] WARNING : No state! Shouldn't happen\n";
                    break;
                }
            }

            // DBGR_INTERMEDIATE_TIME_MEASUREMENT();

            // Tell runtime system how many output items we produced.
            return 0;
        }

        void decodertmp_impl::set_sf(const uint8_t sf) {
            (void) sf;
            std::cerr << "[LoRa Decoder] WARNING : Setting the spreading factor during execution is currently not supported." << std::endl
                      << "Nothing set, kept SF of " << d_sf << "." << std::endl;
        }

        void decodertmp_impl::set_samp_rate(const float samp_rate) {
            (void) samp_rate;
            std::cerr << "[LoRa Decoder] WARNING : Setting the sample rate during execution is currently not supported." << std::endl
                      << "Nothing set, kept SR of " << d_samples_per_second << "." << std::endl;
        }
    } /* namespace lora */
} /* namespace gr */
