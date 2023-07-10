#ifdef ARDUINO_USB_CDC_ON_BOOT
#undef ARDUINO_USB_CDC_ON_BOOT
#endif

// include the library
#include <Arduino.h>
#include <RadioLib.h>
#include <EEPROM.h>

HardwareSerial Serial(0);
// const int EEPROM_SIZE = 1024;
// const int EEPROM_INFO_ADDR = 0;
// struct NodeInfo {
//   int id;
//   uint8_t sf;
//   float bw;
//   float freq;
//   int8_t tp;
//   uint8_t cr;
//   float time_offset;
// };

// SX1278 has the following connections:
// NSS pin:   10
// DIO0 pin:  2
// RESET pin: 9
// DIO1 pin:  3
SX1278 radio = new Module(10, 15, -1, 16);

#if defined(ESP8266) || defined(ESP32)
  ICACHE_RAM_ATTR
#endif

// or using RadioShield
// https://github.com/jgromes/RadioShield
//SX1278 radio = RadioShield.ModuleA;
float fre = 433.5;
float bw = 125.0;
uint8_t control_sf = 7;
uint8_t exp_sf = 10;
uint8_t cr = 5;
uint8_t syncword = 18;
int8_t power = 10;
float symbol_time = pow(2, exp_sf - 7); // uint ms
float esti_symbol = 30;   // estimated symbol number
int delay_max = esti_symbol * symbol_time;
int id;

void setup() {
  Serial.begin(115200);

  // initialize SX1278 with default settings
  Serial.print(F("[SX1278] Initializing ... "));
  int state = radio.begin(fre, bw, exp_sf, cr, syncword, power);
  if (state == RADIOLIB_ERR_NONE) {
    Serial.println(F("success!"));
  } else {
    Serial.print(F("failed, code "));
    Serial.println(state);
    while (true);
  }

  // get id
  // EEPROM.begin(EEPROM_SIZE);
  // NodeInfo info;
  // EEPROM.get(EEPROM_INFO_ADDR, info);
  // id = info.id;
}

void loop() {
  int state = radio.transmit("hello world!testbedmob");
  if (state == RADIOLIB_ERR_NONE) {
    Serial.println(F("transmit success!"));
  }
  delay(5000);

}
