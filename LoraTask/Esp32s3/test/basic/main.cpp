#include <Arduino.h>
#include <EEPROM.h>
#include <Wire.h>
#include <SPI.h>

// IO expander support
#include <PCA9554.h>

// RGB support
#include <FastLED.h>

// LoRa support
#include <RH_RF95.h>

// OLED support
#include <U8g2lib.h>

// Array of leds
CRGB leds[NUM_LEDS];

// LoRa transceiver
RH_RF95 rf95(RA_CS, DIO0);

// IO expander
PCA9554 io(PCA9554A_ADDR);

// OLED
U8G2_SSD1306_128X64_NONAME_F_HW_I2C u8g2(U8G2_R2);

void setup()
{
    // put your setup code here, to run once:

    Serial.begin(115200);
    while (!Serial)
    {
        ;
    }

    Serial.println("src/basic");

    FastLED.addLeds<NEOPIXEL, RGB_BUILTIN>(leds, NUM_LEDS); // GRB ordering is assumed
    FastLED.setBrightness(RGB_BRIGHTNESS / 8);              // Set brightness (0-255)

    if (!rf95.init())
    {
        Serial.println("RF95 init failed");
    }

    Wire.begin();
    io.portMode(ALLOUTPUT);
    io.pinMode(0, INPUT);

    u8g2.begin();
    u8g2.enableUTF8Print();
    u8g2.setFont(u8g2_font_wqy12_t_chinese3);
    // u8g2.setFont(u8g2_font_b12_t_japanese3);
    u8g2.clearBuffer();
    u8g2.setFontDirection(0);
    u8g2.clearBuffer();
    u8g2.setCursor(0, 15);
    u8g2.print("Hello World!");
    u8g2.setCursor(0, 30);
    u8g2.print("你好世界测试程序！"); // Chinese "Hello World"
    //u8g2.print("こんにちは世界"); // Japanese "Hello World"
    u8g2.sendBuffer();
}

void loop()
{
    // put your main code here, to run repeatedly:

    // Serial test
    Serial.print("millis(): ");
    Serial.println(millis());

    // IO expander test
    byte pin = 0;
    io.digitalRead(pin);
    io.digitalWrite(1, pin);
    io.digitalWrite(2, HIGH);
    delay(500);
    io.digitalWrite(2, LOW);
    delay(500);

    // RGB test
    leds[0] = CRGB::Red;
    FastLED.show();
    delay(500);
    leds[0] = CRGB::Black;
    FastLED.show();
    delay(500);

    // LoRa transceiver test
    Serial.println("Sending to rf95_server");
    uint8_t data[] = "Hello World!";
    rf95.send(data, sizeof(data));
    rf95.waitPacketSent();
    uint8_t buf[RH_RF95_MAX_MESSAGE_LEN];
    uint8_t len = sizeof(buf);
    if (rf95.waitAvailableTimeout(2000))
    {
        if (rf95.recv(buf, &len))
        {
            Serial.print("got reply: ");
            Serial.println((char *)buf);
        }
        else
        {
            Serial.println("recv failed");
        }
    }
    else
    {
        Serial.println("No reply, is rf95_server running?");
    }
}
