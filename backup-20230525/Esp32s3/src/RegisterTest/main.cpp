#ifdef ARDUINO_USB_CDC_ON_BOOT
#undef ARDUINO_USB_CDC_ON_BOOT
#endif

#include <Arduino.h>
#include <RadioLib.h>
#include "randomArray.h"

HardwareSerial Serial(0);

// SX1278 has the following connections:
// NSS pin:   10
// DIO0 pin:  2
// RESET pin: 9
// DIO1 pin:  3
SX1278 radio = new Module(10, 15, -1, 16);
// extern unsigned char random_array[1024][1000];

// or using RadioShield
// https://github.com/jgromes/RadioShield
//SX1278 radio = RadioShield.ModuleA;

// flag to indicate that a packet was received
volatile bool transmittedFlag = false;

// flag to indicate frequency must be changed
volatile bool fhssChangeFlag = false;

// the channel frequencies can be generated randomly or hard coded
// NOTE: The frequency list MUST be the same on both sides!
// float channels[] = { 433.5, 433.5 - 0.0625 - 0.009765625 + 0.000244140625, 433.5 + 0.0625 - 0.009765625 + 0.000244140625, 433.775};
// float channels[] = { 433.5, 433.5 - 0.0625, 433.5 + 0.0625, 433.775};
float bw = 0.125;
float gbw = bw/2;
int subchirp_num = 4;
int sf = 10;
int bin = 2<<(sf-1);
int sub_sf = sf - 2 - 2;
int bin_part = bin/subchirp_num;
float C_fre = 433.5;
int downchirp_sync;
int seed_count;
// float channels[] = { 433.5, 433.5 - 0.09375, 433.5 - 0.03125, 433.5 + 0.03125, 433.5 + 0.09375};
float channels[] = { C_fre, C_fre -(bw+gbw)/2-(bw+gbw) - bw/bin, C_fre -(bw+gbw)/2 - bw/bin, C_fre + (bw+gbw)/2 - bw/bin, C_fre + (bw+gbw)/2+(bw+gbw) - bw/bin};
uint8_t channel_jump_array[18*4];
int numberOfChannels = sizeof(channels) / sizeof(float);

// counter to keep track of how many frequency hops were performed
int hopsCompleted = 0;

// counter that increments with each sent packet
int packetCounter = 0;

// save transmission state between loops
int transmissionState = RADIOLIB_ERR_NONE;

int true_bin_count = 0;

// subchirp num


// this is the packet that will be sent
/*
 String longPacket = "Let's create a really long packet to trigger \
 lots of hop interrupts. A packet can be up to 256 bytes long. \
 This packet is 222 bytes so using sf = 9, bw = 125, timeOnAir is \
 1488ms. 1488ms / (9*4.10ms) = 40 hops. Counter: ";
*/
// String longPacket = "1111111";
// uint8_t data[] = {0xFF, 0xFE, 0xFC, 0xF8, 0xF0, 0xE1, 0xC2, 0x85, 0x0B, 0x17, 0x2F, 0x5E, 0xBC, 0x78};
// uint8_t data[] = {0xFF, 0xFE, 0xFC, 0xF8, 0xF0, 0xE1, 0xC2, 0x85, 0x0B, 0x17, 0x2F, 0x5E, 0xBC, 0x78, 0xF1, 0xE3,
//   0xC6, 0x8D, 0x1A, 0x34, 0x68, 0xD0, 0xA0, 0x40, 0x80, 0x01, 0x02, 0x04, 0x08, 0x11, 0x23, 0x47};
uint8_t data[] = {0xFF, 0xFE, 0xFC, 0xF8, 0xF0, 0xE1, 0xC2, 0x85, 0x0B, 0x17, 0x2F, 0x5E, 0xBC, 0x78, 0xF1, 0xE3,
  0xC6, 0x8D, 0x1A, 0x34, 0x68, 0xD0, 0xA0, 0x40, 0x80, 0x01, 0x02, 0x04, 0x08, 0x11, 0x23, 0x47,
  0x8E, 0x1C, 0x38, 0x71, 0xE2, 0xC4, 0x89, 0x12, 0x25, 0x4B, 0x97, 0x2E, 0x5C, 0xB8, 0x70, 0xE0,
  0xC0, 0x81, 0x03, 0x06, 0x0C, 0x19, 0x32, 0x64, 0xC9, 0x92, 0x24, 0x49, 0x93, 0x26, 0x4D, 0x9B,
  0x37, 0x6E, 0xDC, 0xB9, 0x72, 0xE4, 0xC8, 0x90, 0x20, 0x41, 0x82, 0x05, 0x0A, 0x15, 0x2B, 0x56};
int true_bin[] = {0, 56, 112, 168, 224, 280, 336, 392, 448, 504, 560, 616, 672, 728, 784, 840, 896, 952};
// String longPacket = data;

// this function is called when a complete packet
// is transmitted by the module
// IMPORTANT: this function MUST be 'void' type
//            and MUST NOT have any arguments!
#if defined(ESP8266) || defined(ESP32)
  ICACHE_RAM_ATTR
#endif
void setTxFlag(void) {
  transmittedFlag = true;
}

// this function is called when FhssChangeChannel interrupt occurs
// (at the beginning of each transmission)
// IMPORTANT: this function MUST be 'void' type
//            and MUST NOT have any arguments!
#if defined(ESP8266) || defined(ESP32)
  ICACHE_RAM_ATTR
#endif
void setFHSSFlag(void) {
  fhssChangeFlag = true;
}

void setup() {
  Serial.begin(115200);

  // begin radio on home channel
  Serial.print(F("[SX1278] Initializing ... "));
  int state = radio.begin(channels[0], bw*1000, sf, 5, 0x12, 10, 8);
  if (state == RADIOLIB_ERR_NONE) {
    Serial.println(F("success!"));
  } else {
    Serial.print(F("failed, code "));
    Serial.println(state);
    while (true);
  }

  // set hop period in symbols
  // this will also enable FHSS
  state = radio.setFHSSHoppingPeriod(1);
  if (state == RADIOLIB_ERR_NONE) {
    Serial.println(F("success!"));
  } else {
    Serial.print(F("failed, code "));
    Serial.println(state);
    while (true);
  }

  // set the function to call when transmission is finished
  radio.setDio0Action(setTxFlag);

  // set the function to call when we need to hcange frequency
  radio.setDio1Action(setFHSSFlag);

  // set implicit mode and crc disable
  radio.implicitHeader(sizeof(data));
  radio.setCRC(false, false);

  // set random vector
  downchirp_sync = 0;
  // downchirp_sync = random(bin);
  memcpy(channel_jump_array, random_array[downchirp_sync], sizeof(channel_jump_array));
  // randomSeed(downchirp_sync);
  // Serial.print(F("downchirp_sync is: "));
  // Serial.println(downchirp_sync);
  // for(uint32_t i=0; i<sizeof(channel_jump_array) / sizeof(uint8_t); i++){
  //   channel_jump_array[i] = random_array[downchirp_sync][i];
  // }

  Serial.print(channel_jump_array[0]);
  for(uint32_t i=1; i<sizeof(channel_jump_array) / sizeof(uint8_t); i++){
    Serial.print(","); Serial.print(channel_jump_array[i]);
  }
  Serial.println();
  radio.setInvertIQ_upchirp();

  // delay(dev_info.device_id*5*1000); // 根据id间隔5s启动一个
  // start transmitting the first packet
  Serial.print(F("[SX1278] Sending first packet ... "));
  // String packet = longPacket;
  // transmissionState = radio.startTransmit(packet);
  transmissionState = radio.startTransmit(data, sizeof(data), 0U);
}

void loop() {
  // check if the transmission flag is set
  if (transmittedFlag == true) {
    // reset flag
    transmittedFlag = false;

    if (transmissionState == RADIOLIB_ERR_NONE) {
      // packet was successfully sent
      Serial.println(F("transmission finished!"));

    } else {
      Serial.print(F("failed, code "));
      Serial.println(transmissionState);

    }

    // The channel is automatically reset to 0 upon completion
    Serial.print(F("[SX1278] Radio is on channel: "));
    Serial.println(radio.getFHSSChannel());

    // print the number of hops it took
    Serial.print(F("[SX1278] Hops completed: "));
    Serial.println(hopsCompleted);

    // reset the counter
    hopsCompleted = 0;
    true_bin_count = 0;

    // return to home channel before the next transaction
    // radio.setFrequency(channels[radio.getFHSSChannel() % numberOfChannels]);
    radio.setFrequency(channels[0]);
    radio.setSpreadingFactor(sf);
    radio.setBandwidth(bw*1000);
    radio.implicitHeader(sizeof(data));
    radio.setCRC(false, false);
    // add by zkw

    // wait a second before transmitting again
    // delay(300000);
    delay(7000);

    // increment the packet counter
    packetCounter++;

    // send another packet
    Serial.print(F("[SX1278] Sending another packet ... "));
    // randomSeed(1);
    // set random vector
    downchirp_sync = 0;
    // downchirp_sync = random(bin);
    // downchirp_sync = (analogRead(5) % bin) % 1024;
    // random.seed(downchirp_sync);
    memcpy(channel_jump_array, random_array[downchirp_sync], sizeof(channel_jump_array));
    // randomSeed(downchirp_sync);
    // Serial.print(F("downchirp_sync is: "));
    // Serial.println(downchirp_sync);
    // for(uint32_t i=0; i<sizeof(channel_jump_array) / sizeof(uint8_t); i++){
    //   channel_jump_array[i] = random_array[downchirp_sync][i];
    // }

    Serial.print(channel_jump_array[0]);
    for(uint32_t i=1; i<sizeof(channel_jump_array) / sizeof(uint8_t); i++){
      Serial.print(","); Serial.print(channel_jump_array[i]);
    }
    Serial.println();
    radio.setInvertIQ_upchirp();
    transmissionState = radio.startTransmit(data, sizeof(data), 0U);
  }

  // check if we need to do another frequency hop
  if (fhssChangeFlag == true) {
    // we do, change it now
    // int state = radio.setFrequency(channels[radio.getFHSSChannel() % numberOfChannels]);

    // increment the counter
    hopsCompleted++;
    // int state = radio.setFrequency(channels[hopsCompleted % numberOfChannels]);
    // if (state != RADIOLIB_ERR_NONE) {
    //   Serial.print(F("[SX1278] Failed to change frequency, code "));
    //   Serial.println(state);
    // }
    uint32_t jump_part = 3 - true_bin[true_bin_count]/bin_part;
    bool jump_flag = (true_bin[true_bin_count] % bin_part) >= bin_part/2;
    float adjust_bw = bw * (float)true_bin[true_bin_count]/bin;
    // add by zkw
    // add downchirp sync 
    if(hopsCompleted >= 2 && hopsCompleted % 2 == 0){
      // radio.setBandwidthWithoutStandby(bw*1000 / subchirp_num); 
      // radio.setSpreadingFactorWithoutStandby(sub_sf);
      // radio.setInvertIQ_downchirp();
      radio.setInvertIQ_upchirp();
      // jump_flag = downchirp_sync >= bin/2;
      // adjust_bw = bw * (float)downchirp_sync/bin;
      // if(jump_flag)
      //   radio.setFrequency(channels[0] + bw/bin + adjust_bw - bw);
      // else
      //   radio.setFrequency(channels[0] + bw/bin + adjust_bw);
    }
    if(hopsCompleted >= 2 && hopsCompleted % 2 == 1){
      radio.setInvertIQ_upchirp();
    }

    // if(hopsCompleted == 3){
    //   radio.setInvertIQ_upchirp();
    //   radio.setSpreadingFactorWithoutStandby(sub_sf);
    //   if(jump_part == 0 && jump_flag == true)
    //     radio.setFrequency(channels[channel_jump_array[true_bin_count * 4 + 0]] - bw*3/8 - bw + adjust_bw);
    //   else
    //     radio.setFrequency(channels[channel_jump_array[true_bin_count * 4 + 0]] - bw*3/8 + adjust_bw);
      
    // }
    // if(hopsCompleted == 4){
    //   radio.setBandwidthWithoutStandby(bw*1000 / subchirp_num); 
    // }

    // if(hopsCompleted >= 4){
    //   if(hopsCompleted % 4 == 1){
    //     if(jump_part < 2 || (jump_part == 2 && jump_flag == true))
    //       radio.setFrequency(channels[channel_jump_array[true_bin_count * 4 + 2]] + bw/8 - bw + adjust_bw);
    //     else
    //       radio.setFrequency(channels[channel_jump_array[true_bin_count * 4 + 2]] + bw/8 + adjust_bw);
    //   }

    //   else if(hopsCompleted % 4 == 2){
    //     if(jump_part < 3 || (jump_part == 3 && jump_flag == true))
    //       radio.setFrequency(channels[channel_jump_array[true_bin_count * 4 + 3]] + bw*3/8 - bw + adjust_bw);
    //     else
    //       radio.setFrequency(channels[channel_jump_array[true_bin_count * 4 + 3]] + bw*3/8 + adjust_bw);
    //     true_bin_count++;
    //   }

    //   else if(hopsCompleted % 4 == 3){
    //     if(jump_part == 0 && jump_flag == true)
    //       radio.setFrequency(channels[channel_jump_array[true_bin_count * 4 + 0]] - bw*3/8 - bw + adjust_bw); 
    //     else
    //       radio.setFrequency(channels[channel_jump_array[true_bin_count * 4 + 0]] - bw*3/8 + adjust_bw);
    //   }

    //   else if(hopsCompleted % 4 == 0){
    //     if(jump_part < 1 || (jump_part == 1 && jump_flag == true))
    //       radio.setFrequency(channels[channel_jump_array[true_bin_count * 4 + 1]] - bw/8 - bw + adjust_bw);
    //     else
    //       radio.setFrequency(channels[channel_jump_array[true_bin_count * 4 + 1]] - bw/8 + adjust_bw);
    //   }
    // }


    // if(hopsCompleted == 1)
      // radio.setInvertIQ();
    
    // clear the FHSS interrupt
    radio.clearFHSSInt();

    // we're ready to do another hop, clear the flag
    fhssChangeFlag = false;
  }
}
