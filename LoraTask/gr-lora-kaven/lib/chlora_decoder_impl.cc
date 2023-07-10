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
#include "chlora_decoder_impl.h"
#include "tables.h"
#include <string>
#include <pwd.h>
#include <time.h>
#include <cmath>
#include <map>
#include "randomArray.h"
#include <gnuradio/filter/firdes.h>
#include <gnuradio/filter/fir_filter.h>

#include "dbugr.hpp"
#define GRLORA_DEBUG
#define DEBUG false
#define PI acos(-1)
// extern unsigned char random_array[];

namespace gr {
    namespace lora {

        chlora::sptr chlora::make(float samp_rate, uint32_t bandwidth, uint8_t sf, bool implicit, uint8_t cr, bool crc, bool reduced_rate, bool disable_drift_correction, uint32_t pkg_num_all) {
            return gnuradio::get_initial_sptr
                   (new chlora_decoder_impl(samp_rate, bandwidth, sf, implicit, cr, crc, reduced_rate, disable_drift_correction, pkg_num_all));
        }

        /**
         * The private constructor
         */
        chlora_decoder_impl::chlora_decoder_impl(float samp_rate, uint32_t bandwidth, uint8_t sf, bool implicit, uint8_t cr, bool crc, bool reduced_rate, bool disable_drift_correction, uint32_t pkg_num_all)
            : gr::sync_block("chlora_decoder",
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
            d_pkg_num_all = pkg_num_all;
            float leakage_width_array[] = {0.05,0.01,0.015,0.01,0.01,0.01};
            filter_num = 3;
            if(d_sf > 6){
                leakage_width1 = leakage_width_array[d_sf-6];
                leakage_width2 = 1 - leakage_width1;
            }
            preamble_length = 8;
            pos_record.resize(preamble_length);
            if(d_bw == 125e3)
                pass_arg = 0.05;
            else if(d_bw == 250e3)
                pass_arg = 0.075;
            else if(d_bw == 500e3)
                pass_arg = 0.15;
            bw_tmp = d_bw;
            d_gbw = bw_tmp/2;
            carrirFre = {-(bw_tmp+d_gbw)/2-(bw_tmp+d_gbw), -(bw_tmp+d_gbw)/2, (bw_tmp+d_gbw)/2, (bw_tmp+d_gbw)/2+(bw_tmp+d_gbw)};
            build_shift_signal();
            // output_tmp.resize(d_samples_per_symbol);
            // std::cout << carrirFre[0] << " " << carrirFre[1] << " " << carrirFre[2] << " " << carrirFre[3] << std::endl;

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
            d_merge_tmp.resize(d_number_of_bins);
            d_q  = fft_create_plan(d_samples_per_symbol, &d_mult_hf[0], &d_fft[0],     LIQUID_FFT_FORWARD, 0);
            d_qr = fft_create_plan(d_number_of_bins,     &d_tmp[0],     &d_mult_hf[0], LIQUID_FFT_BACKWARD, 0);

            // FFT decoding preparations for zeropadding
            int fftNum_zeropadding = d_number_of_bins * d_decim_factor * pow(2, 10-d_sf);
            int dineNum_zeropadding = d_samples_per_symbol * d_decim_factor * pow(2, 10-d_sf);
            d_fft_zero.resize(dineNum_zeropadding);
            d_mult_hf_zero.resize(dineNum_zeropadding);
            d_tmp_zero.resize(fftNum_zeropadding);
            d_merge_tmp_zero.resize(fftNum_zeropadding);
            d_q_zero  = fft_create_plan(dineNum_zeropadding, &d_mult_hf_zero[0], &d_fft_zero[0],     LIQUID_FFT_FORWARD, 0);
            d_qr_zero = fft_create_plan(dineNum_zeropadding, &d_tmp_zero[0],     &d_mult_hf_zero[0], LIQUID_FFT_BACKWARD, 0);

            // Hamming coding
            fec_scheme fs = LIQUID_FEC_HAMMING84;
            d_h48_fec = fec_create(fs, NULL);

            // Register gnuradio ports
            message_port_register_out(pmt::mp("frames"));
            message_port_register_out(pmt::mp("control"));
            
            // struct  passwd *pwd;
            // pwd = getpwuid(getuid());
            // char a[50] = "rm -rf /mnt/hgfs/share/samples/*";
            // char a[50] = "rm -rf /home/";
            // char b[10] = "/FFT/*";
            // strcat(strcat(a, pwd->pw_name),b);
            // int unused = system(a);
            // unused++;  // supress unused warnning
            timestamp_before = std::time(0);
            build_shift_signal();
            tmp_signal_1.resize(d_samples_per_symbol);
            tmp_signal_2.resize(d_samples_per_symbol);
            tmp_signal_3.resize(d_samples_per_symbol);
            tmp_signal_4.resize(d_samples_per_symbol);
        }

        /**
         * Our virtual destructor.
         */
        chlora_decoder_impl::~chlora_decoder_impl() {
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

        void chlora_decoder_impl::build_ideal_chirps(void) {
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

        void chlora_decoder_impl::values_to_file(const std::string path, const unsigned char *v, const uint32_t length, const uint32_t ppm) {
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

        void chlora_decoder_impl::samples_to_file(const std::string path, const gr_complex *v, const uint32_t length, const uint32_t elem_size) {
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

        void chlora_decoder_impl::samples_to_file_app(const std::string path, const gr_complex *v, const uint32_t length, const uint32_t elem_size) {
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

        void chlora_decoder_impl::samples_debug(const gr_complex *v, const uint32_t length) {
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

        inline void chlora_decoder_impl::instantaneous_frequency(const gr_complex *in_samples, float *out_ifreq, const uint32_t window) {
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

        inline void chlora_decoder_impl::instantaneous_phase(const gr_complex *in_samples, float *out_iphase, const uint32_t window) {
            out_iphase[0] = std::arg(in_samples[0]);

            for (uint32_t i = 1u; i < window; i++) {
                out_iphase[i] = std::arg(in_samples[i]);
                // = the same as atan2(imag(in_samples[i]),real(in_samples[i]));

                // Unwrapped loops from liquid_unwrap_phase
                while ( (out_iphase[i] - out_iphase[i-1]) >  M_PI ) out_iphase[i] -= 2.0f*M_PI;
                while ( (out_iphase[i] - out_iphase[i-1]) < -M_PI ) out_iphase[i] += 2.0f*M_PI;
            }
        }

        float chlora_decoder_impl::cross_correlate_ifreq_fast(const float *samples_ifreq, const float *ideal_chirp, const uint32_t window) {
            float result = 0;
            volk_32f_x2_dot_prod_32f(&result, samples_ifreq, ideal_chirp, window);
            return result;
        }

        float chlora_decoder_impl::cross_correlate_fast(const gr_complex *samples, const gr_complex *ideal_chirp, const uint32_t window) {
            gr_complex result = 0;
            volk_32fc_x2_conjugate_dot_prod_32fc(&result, samples, ideal_chirp, window);
            return abs(result);
        }

        float chlora_decoder_impl::cross_correlate(const gr_complex *samples_1, const gr_complex *samples_2, const uint32_t window) {
            float result = 0.0f;

            for (uint32_t i = 0u; i < window; i++) {
                result += std::real(samples_1[i] * std::conj(samples_2[i]));
            }

            result /= (float)window;

            return result;
        }

        float chlora_decoder_impl::cross_correlate_ifreq(const float *samples_ifreq, const std::vector<float>& ideal_chirp, const uint32_t to_idx) {
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

        void chlora_decoder_impl::fine_sync(const gr_complex* in_samples, int32_t bin_idx, int32_t search_space) {
            int32_t shift_ref = (bin_idx+1) * d_decim_factor;
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
            /*
            if(lag > 0)
                d_fine_sync = std::min(-lag / 2, -1);
            else if(lag < 0)
                d_fine_sync = std::max(-lag / 2, 1);*/

            // Hard limit impact of correction
            /*if(abs(d_fine_sync) >= d_decim_factor / 2)
                d_fine_sync = 0;*/

            //d_fine_sync = 0;
            #ifdef GRLORA_DEBUG
                d_debug << "FINE: " << d_fine_sync << std::endl;
            #endif
        }

        float chlora_decoder_impl::detect_preamble_autocorr(const gr_complex *samples, const uint32_t window) {
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

        float chlora_decoder_impl::determine_energy(const gr_complex *samples) {
            float magsq_chirp[d_samples_per_symbol];
            float energy_chirp = 0;
            volk_32fc_magnitude_squared_32f(magsq_chirp, samples, d_samples_per_symbol);
            volk_32f_accumulator_s32f(&energy_chirp, magsq_chirp, d_samples_per_symbol);

            return energy_chirp;
        }

        void chlora_decoder_impl::determine_snr() {
            if(d_pwr_queue.size() >= 2) {
                float pwr_noise = d_pwr_queue[0];
                float pwr_signal = d_pwr_queue[d_pwr_queue.size()-1];
                d_snr = pwr_signal / pwr_noise;
            }
        }

        float chlora_decoder_impl::detect_downchirp(const gr_complex *samples, const uint32_t window) {
            float samples_ifreq[window];
            instantaneous_frequency(samples, samples_ifreq, window);

            return cross_correlate_ifreq(samples_ifreq, d_downchirp_ifreq, window - 1u);
        }

        float chlora_decoder_impl::detect_upchirp(const gr_complex *samples, const uint32_t window, int32_t *index) {
            float samples_ifreq[window*2];
            instantaneous_frequency(samples, samples_ifreq, window*2);

            return sliding_norm_cross_correlate_upchirp(samples_ifreq, window, index);
        }

        float chlora_decoder_impl::sliding_norm_cross_correlate_upchirp(const float *samples_ifreq, const uint32_t window, int32_t *index) {
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

        float chlora_decoder_impl::stddev(const float *values, const uint32_t len, const float mean) {
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
        uint32_t chlora_decoder_impl::get_shift_fft(const gr_complex *samples) {
            float fft_mag[d_number_of_bins];

            samples_to_file("/tmp/data", &samples[0], d_samples_per_symbol, sizeof(gr_complex));

            // Multiply with ideal downchirp
            for (uint32_t i = 0u; i < d_samples_per_symbol; i++) {
                d_mult_hf[i] = samples[i] * d_downchirp[i];
            }

            samples_to_file("/tmp/mult", &d_mult_hf[0], d_samples_per_symbol, sizeof(gr_complex));

            // Perform FFT
            fft_execute(d_q);

            // Decimate. Note: assumes fft size is multiple of decimation factor and number of bins is even
            // This decimation should be identical to numpy's approach
            const uint32_t N = d_number_of_bins;
            memcpy(&d_tmp[0],               &d_fft[0],                                     (N + 1u) / 2u * sizeof(gr_complex));
            memcpy(&d_tmp[ (N + 1u) / 2u ], &d_fft[d_samples_per_symbol - (N / 2u)],        N / 2u * sizeof(gr_complex));
            d_tmp[N / 2u] += d_fft[N / 2u];

            // Get magnitude
            for (uint32_t i = 0u; i < d_number_of_bins; i++) {
                fft_mag[i] = std::abs(d_tmp[i]);
            }

            samples_to_file("/tmp/fft", &d_tmp[0], d_number_of_bins, sizeof(gr_complex));

            fft_execute(d_qr); // For debugging
            samples_to_file("/tmp/resampled", &d_mult_hf[0], d_number_of_bins, sizeof(gr_complex));

            // Return argmax here
            return (std::max_element(fft_mag, fft_mag + d_number_of_bins) - fft_mag);
        }

        uint32_t chlora_decoder_impl::max_frequency_gradient_idx(const gr_complex *samples) {
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

        bool chlora_decoder_impl::demodulate(const gr_complex *samples, const bool is_first) {
            // DBGR_TIME_MEASUREMENT_TO_FILE("SFxx_method");
            bool reduced_rate = is_first || d_reduced_rate;

            // DBGR_START_TIME_MEASUREMENT(false, "only");

            uint32_t bin_idx = max_frequency_gradient_idx(samples);
            //uint32_t bin_idx = get_shift_fft(samples);
            if(d_enable_fine_sync)
                fine_sync(samples, bin_idx, std::max(d_decim_factor / 4u, 2u));

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
        void chlora_decoder_impl::deinterleave(const uint32_t ppm) {
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

        void chlora_decoder_impl::decode(const bool is_header) {
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

        void chlora_decoder_impl::msg_lora_frame(void) {
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

        void chlora_decoder_impl::deshuffle(const uint8_t *shuffle_pattern, const bool is_header) {
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

        void chlora_decoder_impl::dewhiten(const uint8_t *prng) {
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

        void chlora_decoder_impl::hamming_decode(bool is_header) {
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
        void chlora_decoder_impl::hamming_decode_soft(bool is_header) {
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

        void chlora_decoder_impl::extract_data_only(bool is_header) {
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
        void chlora_decoder_impl::determine_cfo(const gr_complex *samples) {
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
        float chlora_decoder_impl::experimental_determine_cfo(const gr_complex *samples, uint32_t window) {
            gr_complex mult[window];
            float mult_ifreq[window];

            volk_32fc_x2_multiply_32fc(mult, samples, &d_downchirp[0], window);
            instantaneous_frequency(mult, mult_ifreq, window);

            return mult_ifreq[256] / (2.0 * M_PI) * d_samples_per_second;
        }

        float chlora_decoder_impl::detect_amp_rate(const gr_complex *samples, const gr_complex *ideal_chirp) {
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

        bool chlora_decoder_impl::detect_preamble(const gr_complex *samples, uint32_t condition) {
            // save fft_merge_abs value
            float fft_mag[d_number_of_bins];
            // dechirp
            for (uint32_t i = 0u; i < d_samples_per_symbol; i++) {
                d_mult_hf[i] = samples[i] * d_downchirp[i];
            }
            // Perform FFT
            fft_execute(d_q);
            // merge fft
            memcpy(&d_tmp[0], &d_fft[0], d_number_of_bins * sizeof(gr_complex));
            memcpy(&d_merge_tmp[0], &d_fft[d_samples_per_symbol-d_number_of_bins], d_number_of_bins * sizeof(gr_complex));
            for (uint32_t i = 0u; i < d_number_of_bins; i++) {
                d_tmp[i] = d_tmp[i] + d_merge_tmp[i];
                // get amptitude
                fft_mag[i] = std::abs(d_tmp[i]);
            }
            // find max amptitude
            float max_amptitude = *std::max_element(fft_mag, fft_mag + d_number_of_bins);
            // compute mean
            float sum_amptitude = 0.0;
            for (uint32_t i = 0u; i < d_number_of_bins; i++) {
                sum_amptitude += fft_mag[i]; 
            }
            float mean_amptitude = sum_amptitude / d_number_of_bins;
            // std::cout << max_amptitude << " " << mean_amptitude << std::endl; 
            if(max_amptitude > condition*mean_amptitude)
                return true;
            else return false;
        }

        std::vector<int> chlora_decoder_impl::get_filter_sorted_peak(const gr_complex *samples) {
            // save fft_merge_abs value
            int fftNum_zeropadding = d_number_of_bins * d_decim_factor * pow(2, 10-d_sf);
            int dineNum_zeropadding = d_samples_per_symbol * d_decim_factor * pow(2, 10-d_sf);
            std::vector<float> fft_mag;
            fft_mag.resize(fftNum_zeropadding);
            // float fft_mag[fftNum_zeropadding];
            // dechirp
            if(DEBUG)
                std::cout << "dechirp before" << std::endl;
            for (uint32_t i = 0u; i < d_samples_per_symbol; i++) {
                d_mult_hf_zero[i] = samples[i] * d_downchirp[i];
            }
            // add zero
            for (uint32_t i = d_samples_per_symbol; i < dineNum_zeropadding; i++) {
                d_mult_hf_zero[i] = 0;
            }
            // Perform FFT
            if(DEBUG)
                std::cout << "fft before" << std::endl;
            fft_execute(d_q_zero);
            if(DEBUG)
                std::cout << "fft after" << std::endl;
            // merge fft
            memcpy(&d_tmp_zero[0], &d_fft_zero[0], fftNum_zeropadding * sizeof(gr_complex));
            memcpy(&d_merge_tmp_zero[0], &d_fft_zero[dineNum_zeropadding - fftNum_zeropadding], fftNum_zeropadding * sizeof(gr_complex));
            for (uint32_t i = 0u; i < fftNum_zeropadding; i++) {
                d_tmp_zero[i] = d_tmp_zero[i] + d_merge_tmp_zero[i];
                // get amptitude
                fft_mag[i] = std::abs(d_tmp_zero[i]);
            }
            if(DEBUG)
                std::cout << "get amptitude" << std::endl;
            // sort amptitude index
            std::vector<int> index(fftNum_zeropadding);
            for (int i = 0 ; i != fftNum_zeropadding ; i++) {
                index[i] = i;
            }
            sort(index.begin(), index.end(), [&](const int& a, const int& b) {
                return (fft_mag[a] > fft_mag[b]);
            });
            if(DEBUG)
                std::cout << "sort!" << std::endl;
            // pick the filter_num peak which sidelobe is eliminated
            std::vector<int> pos(3);
            pos[0] = index[0];
            int count = 1;
            for(int i=1; i<fftNum_zeropadding; i++){
                int temp1 = index[i];
                for(int k=0; k<count; k++){
                    int temp2 = pos[k];
                    float dif = abs(temp2 - temp1);
                    if(dif < fftNum_zeropadding*leakage_width1 || dif > fftNum_zeropadding*leakage_width2) // the same as sidelobe
                        break;
                    if(k == count - 1){ // pick the order peak
                        pos[count] = temp1;
                        count++;
                        break;
                    }
                }
                if(count == filter_num)
                    break;
            }
            if(DEBUG)
                std::cout << "find peak" << std::endl;
            return pos;
        }

        int chlora_decoder_impl::get_SFD_bin(const gr_complex *samples){
            // save fft_merge_abs value
            int fftNum_zeropadding = d_number_of_bins * d_decim_factor * pow(2, 10-d_sf);
            int dineNum_zeropadding = d_samples_per_symbol * d_decim_factor * pow(2, 10-d_sf);
            float fft_mag[fftNum_zeropadding];
            // dechirp
            for (uint32_t i = 0u; i < d_samples_per_symbol; i++) {
                d_mult_hf_zero[i] = samples[i] * d_upchirp[i];
            }
            // add zero
            for (uint32_t i = d_samples_per_symbol; i < dineNum_zeropadding; i++) {
                d_mult_hf_zero[i] = 0;
            }
            // Perform FFT
            fft_execute(d_q_zero);
            // merge fft
            memcpy(&d_tmp_zero[0], &d_fft_zero[0], fftNum_zeropadding * sizeof(gr_complex));
            memcpy(&d_merge_tmp_zero[0], &d_fft_zero[dineNum_zeropadding - fftNum_zeropadding], fftNum_zeropadding * sizeof(gr_complex));
            for (uint32_t i = 0u; i < fftNum_zeropadding; i++) {
                d_tmp_zero[i] = d_tmp_zero[i] + d_merge_tmp_zero[i];
                // get amptitude
                fft_mag[i] = std::abs(d_tmp_zero[i]);
            }
            return std::max_element(fft_mag, fft_mag + fftNum_zeropadding) - fft_mag;
        }

        int chlora_decoder_impl::get_preamble_bin(std::vector<std::vector<int>> array){
            std::vector<int> candidate(array.size()*array[0].size());
            int fftNum_zeropadding = d_number_of_bins * d_decim_factor * pow(2, 10-d_sf);
            int upchirp_ref = 0;
            if(array[0][0] == array[1][0])
                upchirp_ref = array[0][0];
            else upchirp_ref = array[1][0]; 
            // std::cout << "upchirp_ref: " << upchirp_ref << std::endl;
            int count = 0;
            for(int i=0; i<array.size(); i++){
                for(int k=0; k<array[0].size(); k++){
                    int dif = std::abs(array[i][k] - upchirp_ref);
                    if(dif < fftNum_zeropadding*leakage_width1 || dif > fftNum_zeropadding*leakage_width2){
                        candidate[count] = array[i][k];
                        count++;
                    }
                }
            }
            // for(int i=0; i<count; i++){
            //     std::cout << candidate[i] << " ";
            // }
            candidate.resize(count);
            // create map to count the most frequency number
            std::map<int, int> map;
            int max_fre_num = 0;
            for(int i=0; i<count; i++){
                map[candidate[i]]++;
                if(map[candidate[i]] > map[max_fre_num])
                    max_fre_num = candidate[i];
            }
            int count_preamble_num = 0;
            for(int i=0; i<count; i++){
                int dif = std::abs(candidate[i] - max_fre_num);
                if( dif < 4 || dif > fftNum_zeropadding-4)
                    count_preamble_num++;
            }
            premable_count = count_preamble_num;
            return max_fre_num;
        }

        void chlora_decoder_impl::caculate_cfo_winoff(int upchirp_bin, int downchirp_bin){
            float upbin = upchirp_bin, downbin = downchirp_bin;
            float fftNum_zeropadding = d_number_of_bins * d_decim_factor * pow(2, 10-d_sf);
            float cfo_bin = 0;
            if(upbin + downbin < fftNum_zeropadding*0.5){
                cfo_bin = upbin + downbin - 2;
                cfo = -cfo_bin/2/fftNum_zeropadding * d_bw;
                windows_offset = (downbin - upbin) / pow(2, 11-d_sf);
            }else if(upbin + downbin > fftNum_zeropadding*1.5){
                cfo_bin = upbin + downbin - fftNum_zeropadding*2*((float)d_decim_factor/8) - 2;
                cfo = -cfo_bin/2/fftNum_zeropadding * d_bw;
                windows_offset = (downbin - upbin) / pow(2, 11-d_sf);
            }else{
                cfo_bin = upbin + downbin - fftNum_zeropadding - 2;
                cfo = -cfo_bin/2/fftNum_zeropadding * d_bw;
                windows_offset = (fftNum_zeropadding - (upbin - downbin)) / pow(2, 11-d_sf);
            }
        }

        void chlora_decoder_impl::test_method(const gr_complex *samples){
            // save fft_merge_abs value
            int fftNum_zeropadding = d_number_of_bins * d_decim_factor * pow(2, 10-d_sf);
            int dineNum_zeropadding = d_samples_per_symbol * d_decim_factor * pow(2, 10-d_sf);
            float fft_mag[fftNum_zeropadding];
            // dechirp
            for (uint32_t i = 0u; i < d_samples_per_symbol; i++) {
                d_mult_hf_zero[i] = samples[i] * d_downchirp[i];
            }
            // add zero
            for (uint32_t i = d_samples_per_symbol; i < dineNum_zeropadding; i++) {
                d_mult_hf_zero[i] = 0;
            }
            // Perform FFT
            fft_execute(d_q_zero);
            // merge fft
            memcpy(&d_tmp_zero[0], &d_fft_zero[0], fftNum_zeropadding * sizeof(gr_complex));
            memcpy(&d_merge_tmp_zero[0], &d_fft_zero[dineNum_zeropadding - fftNum_zeropadding], fftNum_zeropadding * sizeof(gr_complex));
            for (uint32_t i = 0u; i < fftNum_zeropadding; i++) {
                d_tmp_zero[i] = d_tmp_zero[i] + d_merge_tmp_zero[i];
                // get amptitude
                fft_mag[i] = std::abs(d_tmp_zero[i]);
            }
            std::cout << std::max_element(fft_mag, fft_mag + fftNum_zeropadding) - fft_mag << std::endl;
        }

        void chlora_decoder_impl::build_cfo_idealchirp(void){
            d_downchirp_cfo.resize(d_samples_per_symbol);
            d_upchirp_cfo.resize(d_samples_per_symbol);
            
            const double T       = -0.5 * d_bw * d_symbols_per_second;
            const double f0      = (d_bw / 2.0) + cfo;
            const double f1      = (d_bw / 2.0) - cfo;
            const double pre_dir = 2.0 * M_PI;
            double t;
            gr_complex cmx       = gr_complex(1.0f, 1.0f);

            for (uint32_t i = 0u; i < d_samples_per_symbol; i++) {
                // Width in number of samples = samples_per_symbol
                // See https://en.wikipedia.org/wiki/Chirp#Linear
                t = d_dt * i;
                d_downchirp_cfo[i] = cmx * gr_expj(pre_dir * t * (f0 + T * t));
                d_upchirp_cfo[i]   = cmx * gr_expj(pre_dir * t * (f1 + T * t) * -1.0f);
            }
        }

        int chlora_decoder_impl::get_upchirp_bin_cfo(const gr_complex *samples) {
            // save fft_merge_abs value
            float fft_mag[d_number_of_bins];
            // dechirp
            for (uint32_t i = 0u; i < d_samples_per_symbol; i++) {
                d_mult_hf[i] = samples[i] * d_downchirp_cfo[i];
            }
            // Perform FFT
            fft_execute(d_q);
            // merge fft
            memcpy(&d_tmp[0], &d_fft[0], d_number_of_bins * sizeof(gr_complex));
            memcpy(&d_merge_tmp[0], &d_fft[d_samples_per_symbol-d_number_of_bins], d_number_of_bins * sizeof(gr_complex));
            for (uint32_t i = 0u; i < d_number_of_bins; i++) {
                d_tmp[i] = d_tmp[i] + d_merge_tmp[i];
                // get amptitude
                fft_mag[i] = std::abs(d_tmp[i]);
            }
            // find max amptitude
            return std::max_element(fft_mag, fft_mag + d_number_of_bins) - fft_mag;
        }

        int chlora_decoder_impl::get_downchirp_bin_cfo(const gr_complex *samples) {
            // save fft_merge_abs value
            float fft_mag[d_number_of_bins];
            // dechirp
            for (uint32_t i = 0u; i < d_samples_per_symbol; i++) {
                d_mult_hf[i] = samples[i] * d_upchirp_cfo[i];
            }
            // Perform FFT
            fft_execute(d_q);
            // merge fft
            memcpy(&d_tmp[0], &d_fft[0], d_number_of_bins * sizeof(gr_complex));
            memcpy(&d_merge_tmp[0], &d_fft[d_samples_per_symbol-d_number_of_bins], d_number_of_bins * sizeof(gr_complex));
            for (uint32_t i = 0u; i < d_number_of_bins; i++) {
                d_tmp[i] = d_tmp[i] + d_merge_tmp[i];
                // get amptitude
                fft_mag[i] = std::abs(d_tmp[i]);
            }
            // find max amptitude
            return std::max_element(fft_mag, fft_mag + d_number_of_bins) - fft_mag;
        }

        void chlora_decoder_impl::create_jumpchannel_array(unsigned int num, unsigned int channel_num){
            u_int64_t temp = num;
            int temp_array[] = {11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61};
            for(int k=0; k<10; k++){
                u_int64_t temp1 = (temp << 12) * temp_array[k];
                for(int i=0; i<10; i++){
                    channel_jump_array[k*10+i] = temp1 % channel_num;
                    temp1 = temp1 / 4;
                }
            }
        }
        
        void chlora_decoder_impl::build_shift_signal(){
            shift_cos_1.resize(d_samples_per_symbol);
            shift_cos_2.resize(d_samples_per_symbol);
            shift_cos_3.resize(d_samples_per_symbol);
            shift_cos_4.resize(d_samples_per_symbol);
            const double pre_dir = 2.0 * M_PI;
            double t;

            for (uint32_t i = 0u; i < d_samples_per_symbol; i++) {
                t = d_dt * i;
                shift_cos_1[i] = gr_expj(pre_dir * t * -carrirFre[0]);
                shift_cos_2[i] = gr_expj(pre_dir * t * -carrirFre[1]);
                shift_cos_3[i] = gr_expj(pre_dir * t * -carrirFre[2]);
                shift_cos_4[i] = gr_expj(pre_dir * t * -carrirFre[3]);
            }
        }

        void chlora_decoder_impl::shift_signal(const gr_complex *samples){
            for (uint32_t i = 0u; i < d_samples_per_symbol; i++) {
                tmp_signal_1[i] = samples[i] * shift_cos_1[i] * 2.0F;
                tmp_signal_2[i] = samples[i] * shift_cos_2[i] * 2.0F;
                tmp_signal_3[i] = samples[i] * shift_cos_3[i] * 2.0F;
                tmp_signal_4[i] = samples[i] * shift_cos_4[i] * 2.0F;
            }
            // std::cout << "start" << std::endl;
            gr_complex output_tmp[d_samples_per_symbol];
            apply_lowpass_filter(&tmp_signal_1[0], output_tmp);
            memcpy(&tmp_signal_1[0], &output_tmp[0], d_samples_per_symbol * sizeof(gr_complex));
            apply_lowpass_filter(&tmp_signal_2[0], output_tmp);
            memcpy(&tmp_signal_2[0], &output_tmp[0], d_samples_per_symbol * sizeof(gr_complex));
            apply_lowpass_filter(&tmp_signal_3[0], output_tmp);
            memcpy(&tmp_signal_3[0], &output_tmp[0], d_samples_per_symbol * sizeof(gr_complex));
            apply_lowpass_filter(&tmp_signal_4[0], output_tmp);
            memcpy(&tmp_signal_4[0], &output_tmp[0], d_samples_per_symbol * sizeof(gr_complex));
        }

        int chlora_decoder_impl::align_windows(){
            int record[4];
            for(int subchirpCount=0; subchirpCount<4; subchirpCount++){
                // save fft_merge_abs value
                int fftNum_zeropadding = d_number_of_bins * d_decim_factor * pow(2, 10-d_sf);
                int dineNum_zeropadding = d_samples_per_symbol * d_decim_factor * pow(2, 10-d_sf);
                float fft_mag[d_number_of_bins];
                // dechirp
                std::vector<gr_complex> samples_tmp;
                // std::cout << (int)random_array[upchirp_random][subchirpCount] << std::endl;
                switch(random_array[upchirp_random][subchirpCount]){
                    case 1:  samples_tmp = tmp_signal_1; break;
                    case 2:  samples_tmp = tmp_signal_2; break;
                    case 3:  samples_tmp = tmp_signal_3; break;
                    case 4:  samples_tmp = tmp_signal_4; break;
                }
                std::fill(d_mult_hf_zero.begin(), d_mult_hf_zero.end(), 0);
                for (uint32_t i = 0; i < d_samples_per_symbol/4; i++) {
                    d_mult_hf_zero[i] = samples_tmp[i+subchirpCount*d_samples_per_symbol/4] * d_downchirp_cfo[i+subchirpCount*d_samples_per_symbol/4];
                }
                // Perform FFT
                fft_execute(d_q_zero);
                // merge fft
                memcpy(&d_tmp_zero[0], &d_fft_zero[0], d_number_of_bins * sizeof(gr_complex));
                memcpy(&d_merge_tmp_zero[0], &d_fft_zero[d_samples_per_symbol - d_number_of_bins], d_number_of_bins * sizeof(gr_complex));
                for (uint32_t i = 0u; i < d_number_of_bins; i++) {
                    fft_mag[i] = std::abs(d_tmp_zero[i] + d_merge_tmp_zero[i]);
                }
                record[subchirpCount] = std::max_element(fft_mag, fft_mag + d_number_of_bins) - fft_mag;
                // std::cout << record[subchirpCount] << " ";
                // std::cout << std::endl;
            }
            // create map to count the most frequency number
            std::map<int, int> map;
            int max_fre_num = 0;
            for(int i=0; i<4; i++){
                map[record[i]]++;
                if(map[record[i]] > map[max_fre_num])
                    max_fre_num = record[i];
            }
            std::cout << "align windows: " << max_fre_num << std::endl;
            return max_fre_num;
        }

        void chlora_decoder_impl::demodulate(){
            for(int subchirpCount=0; subchirpCount<4; subchirpCount++){
                // save fft_merge_abs value
                float fft_mag[d_number_of_bins];
                // dechirp
                std::vector<gr_complex> samples_tmp;
                // std::cout << (int)random_array[upchirp_random][subchirpCount] << std::endl;
                switch(random_array[upchirp_random][pos_record_count*4 + subchirpCount]){
                    case 1:  samples_tmp = tmp_signal_1; break;
                    case 2:  samples_tmp = tmp_signal_2; break;
                    case 3:  samples_tmp = tmp_signal_3; break;
                    case 4:  samples_tmp = tmp_signal_4; break;
                }
                std::fill(d_mult_hf.begin(), d_mult_hf.end(), 0);
                for (uint32_t i = 0; i < d_samples_per_symbol/4; i++) {
                    d_mult_hf[i] = samples_tmp[i+subchirpCount*d_samples_per_symbol/4] * d_downchirp_cfo[i+subchirpCount*d_samples_per_symbol/4];
                }
                // Perform FFT
                fft_execute(d_q);
                // merge fft
                memcpy(&d_tmp[0], &d_fft[0], d_number_of_bins * sizeof(gr_complex));
                memcpy(&d_merge_tmp[0], &d_fft[d_samples_per_symbol - d_number_of_bins], d_number_of_bins * sizeof(gr_complex));
                for (uint32_t i = 0u; i < d_number_of_bins; i++) {
                    fft_mag[i] = std::abs(d_tmp[i] + d_merge_tmp[i]);
                }
                int tmp = std::max_element(fft_mag, fft_mag + d_number_of_bins) - fft_mag;
                std::cout << tmp << " ";
            }
            std::cout << std::endl;
        }

        void chlora_decoder_impl::apply_lowpass_filter(const gr_complex *samples, gr_complex *output)
        {
            // Compute filter taps
            const std::vector<float> taps = gr::filter::firdes::low_pass(1, d_samples_per_second, d_samples_per_second*0.05, d_samples_per_second*0.2);
            // std::cout << "pahse 1" << std::endl;
            // Create filter object and apply it to the input vector
            gr::filter::kernel::fir_filter_ccf filter(1, taps);
            // gr_complex output[d_samples_per_symbol];
            // std::cout << "pahse 2" << std::endl;
            filter.filterN(output, samples, (unsigned long)d_samples_per_symbol);
            // std::cout << "pahse 3" << std::endl;
        }

        int chlora_decoder_impl::detect_active_channel(){
            int max_amp = 0;
            int active_channel = 0;
            for(int channel=0; channel<4; channel++){
                // save fft_merge_abs value
                float fft_mag[d_number_of_bins];
                // dechirp
                std::vector<gr_complex> samples_tmp;
                // std::cout << (int)random_array[upchirp_random][subchirpCount] << std::endl;
                // std::cout << "pahse 1" << std::endl;
                switch(channel){
                    case 0:  samples_tmp = tmp_signal_1; break;
                    case 1:  samples_tmp = tmp_signal_2; break;
                    case 2:  samples_tmp = tmp_signal_3; break;
                    case 3:  samples_tmp = tmp_signal_4; break;
                }
                // memcpy(&d_tmp[0], &d_fft[0], d_number_of_bins * sizeof(gr_complex));
                // std::fill(d_mult_hf.begin(), d_mult_hf.end(), 0);
                // std::cout << "pahse 2" << std::endl;
                for (uint32_t i = 0; i < d_samples_per_symbol; i++) {
                    d_mult_hf[i] = samples_tmp[i] * d_downchirp[i];
                }
                // Perform FFT
                fft_execute(d_q);
                // std::cout << "pahse 3" << std::endl;
                // merge fft
                memcpy(&d_tmp[0], &d_fft[0], d_number_of_bins * sizeof(gr_complex));
                memcpy(&d_merge_tmp[0], &d_fft[d_samples_per_symbol - d_number_of_bins], d_number_of_bins * sizeof(gr_complex));
                for (uint32_t i = 0u; i < d_number_of_bins; i++) {
                    fft_mag[i] = std::abs(d_tmp[i] + d_merge_tmp[i]);
                }
                // std::cout << "pahse 4" << std::endl;
                float amp = *std::max_element(fft_mag, fft_mag + d_number_of_bins);
                if(amp > max_amp){
                    max_amp = amp;
                    active_channel = channel;
                }
                // std::cout << "max_amp: " << max_amp << std::endl;
            }
            return active_channel;
        }

        int chlora_decoder_impl::work(int noutput_items,
                               gr_vector_const_void_star& input_items,
                               gr_vector_void_star&       output_items) {
            (void) noutput_items;
            (void) output_items;
            
            struct  passwd *pwd;
            pwd = getpwuid(getuid());

            const gr_complex *input     = (gr_complex *) input_items[0];
            //const gr_complex *raw_input = (gr_complex *) input_items[1]; // Input bypassed by low pass filter

            d_fine_sync = 0; // Always reset fine sync
            // uint32_t timestamp_after = std::time(0);
            // std::cout << timestamp_after - timestamp_before << std::endl;
            // if(timestamp_after - timestamp_before > 60)
            //     exit(0);

            switch (d_state) {
                case gr::lora::DecoderState::DETECT: {
                    float correlation = detect_preamble_autocorr(input, d_samples_per_symbol);
                    if(correlation >= 0.90f) {
                        shift_signal(input + d_samples_per_symbol);
                        // std::cout << "pahse 1" << std::endl;
                        preamble_channel = detect_active_channel();
                        std::cout << "Detect pass! active channel: " << preamble_channel << std::endl;
                        d_state = gr::lora::DecoderState::SYNC;
                        pos_record_count = 0;
                        break;
                    }
                    consume_each(d_samples_per_symbol);
                    break;
                }

                case gr::lora::DecoderState::SYNC: {
                    if(pos_record_count < preamble_length){
                        shift_signal(input);
                        // std::cout << "pahse 1" << std::endl;
                        gr_complex signal_tmp[d_samples_per_symbol];
                        switch(preamble_channel){
                            case 0:  memcpy(&signal_tmp[0], &tmp_signal_1[0], d_samples_per_symbol * sizeof(gr_complex)); break;
                            case 1:  memcpy(&signal_tmp[0], &tmp_signal_2[0], d_samples_per_symbol * sizeof(gr_complex)); break;
                            case 2:  memcpy(&signal_tmp[0], &tmp_signal_3[0], d_samples_per_symbol * sizeof(gr_complex)); break;
                            case 3:  memcpy(&signal_tmp[0], &tmp_signal_4[0], d_samples_per_symbol * sizeof(gr_complex)); break;
                        }
                        // std::cout << "pahse 2" << std::endl;
                        pos_record[pos_record_count] = get_filter_sorted_peak(signal_tmp);
                        // get_filter_sorted_peak(input);
                        // detect_preamble(input, 3);
                        // test_method(input);
                        // int downchip_bin = get_SFD_bin(input);
                        // std::cout << "downchirp_bin: " << downchip_bin << std::endl;
                    }else if(pos_record_count == preamble_length + 3){
                        shift_signal(input);
                        gr_complex signal_tmp[d_samples_per_symbol];
                        switch(preamble_channel){
                            case 0:  memcpy(&signal_tmp[0], &tmp_signal_1[0], d_samples_per_symbol * sizeof(gr_complex)); break;
                            case 1:  memcpy(&signal_tmp[0], &tmp_signal_2[0], d_samples_per_symbol * sizeof(gr_complex)); break;
                            case 2:  memcpy(&signal_tmp[0], &tmp_signal_3[0], d_samples_per_symbol * sizeof(gr_complex)); break;
                            case 3:  memcpy(&signal_tmp[0], &tmp_signal_4[0], d_samples_per_symbol * sizeof(gr_complex)); break;
                        }
                        int downchip_bin = get_SFD_bin(signal_tmp);
                        std::cout << "downchirp_bin: " << downchip_bin << std::endl;
                        for(int i=0; i<pos_record.size(); i++){
                            std::vector<int> temp = pos_record[i];
                            for(int k=0; k<temp.size(); k++){
                                std::cout << temp[k] << " ";
                            }
                            std::cout << std::endl;
                        }
                        int upchirp_bin = get_preamble_bin(pos_record);
                        std::cout << "upchirp_bin: " << upchirp_bin << std::endl;
                        std::cout << "preamble_count: " << premable_count << std::endl;
                        caculate_cfo_winoff(upchirp_bin, downchip_bin);
                        std::cout << "cfo: " << cfo << "   " << "winoff: " << windows_offset << std::endl;
                        build_cfo_idealchirp();
                        pos_record_count = 0;
                        consume_each(windows_offset + 0.25 * d_samples_per_symbol);
                        // if(windows_offset > 0)
                        //     consume_each(d_samples_per_symbol - windows_offset);
                        // else consume_each(d_samples_per_symbol - windows_offset);
                        d_state = gr::lora::DecoderState::FIND_SFD;
                        break;
                    }
                    pos_record_count++;
                    consume_each(d_samples_per_symbol);
                    break;
                }

                case gr::lora::DecoderState::FIND_SFD: {
                    if (pos_record_count == 2) { 
                        shift_signal(input);
                        gr_complex signal_tmp[d_samples_per_symbol];
                        switch(preamble_channel){
                            case 0:  memcpy(&signal_tmp[0], &tmp_signal_1[0], d_samples_per_symbol * sizeof(gr_complex)); break;
                            case 1:  memcpy(&signal_tmp[0], &tmp_signal_2[0], d_samples_per_symbol * sizeof(gr_complex)); break;
                            case 2:  memcpy(&signal_tmp[0], &tmp_signal_3[0], d_samples_per_symbol * sizeof(gr_complex)); break;
                            case 3:  memcpy(&signal_tmp[0], &tmp_signal_4[0], d_samples_per_symbol * sizeof(gr_complex)); break;
                        }
                        int upchirp_bin = get_downchirp_bin_cfo(signal_tmp);
                        std::cout << "upchirp_bin: " << upchirp_bin << std::endl;
                        // memcpy(channel_jump_array, random_array[upchirp_bin], sizeof(channel_jump_array));
                        // create_jumpchannel_array(upchirp_bin, 4);
                        // for(int i=0; i<100; i++){
                        //     std::cout << (int)random_array[upchirp_bin][i] << " ";
                        // }
                        // std::cout << std::endl;
                        upchirp_random = upchirp_bin;
                        // consume_each(d_samples_per_symbol);
                        d_state = gr::lora::DecoderState::SYNC_SUB;
                        break;
                    }
                    pos_record_count++;
                    consume_each(d_samples_per_symbol);
                    break;
                }

                case gr::lora::DecoderState::SYNC_SUB: {
                    
                    shift_signal(input + d_samples_per_symbol);
                    int tmp = align_windows();
                    consume_each(d_samples_per_symbol - tmp);
                    d_state = gr::lora::DecoderState::DECODE;
                    pos_record_count = 0;
                    break;
                }

                case gr::lora::DecoderState::DECODE: {
                    if(pos_record_count < 16){
                        shift_signal(input);
                        demodulate();
                    }else exit(0);
                    pos_record_count++;
                    consume_each(d_samples_per_symbol);
                    break;
                }

                default: {
                    std::cerr << "[LoRa CLIPER] WARNING : No state! Shouldn't happen\n";
                    break;
                }
            }

            // DBGR_INTERMEDIATE_TIME_MEASUREMENT();

            // Tell runtime system how many output items we produced.
            return 0;
        }

        void chlora_decoder_impl::set_sf(const uint8_t sf) {
            (void) sf;
            std::cerr << "[LoRa Decoder] WARNING : Setting the spreading factor during execution is currently not supported." << std::endl
                      << "Nothing set, kept SF of " << d_sf << "." << std::endl;
        }

        void chlora_decoder_impl::set_samp_rate(const float samp_rate) {
            (void) samp_rate;
            std::cerr << "[LoRa Decoder] WARNING : Setting the sample rate during execution is currently not supported." << std::endl
                      << "Nothing set, kept SR of " << d_samples_per_second << "." << std::endl;
        }
    } /* namespace lora */
} /* namespace gr */
