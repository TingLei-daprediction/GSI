#include <stdio.h>
#include <string.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#define POLYNOMIAL 0xEDB88320

//# Static lookup table initialized at compile-time
 static uint32_t crc32_table[256] = {0};

 // Function to initialize the CRC32 lookup table (if not already initialized)
 static void init_crc32_table() {
     if (crc32_table[1] != 0) return;  // Prevent redundant initialization

         for (uint32_t i = 0; i < 256; i++) {
                 uint32_t crc = i;
                 for (int j = 0; j < 8; j++) {
                     if (crc & 1)
                       crc = (crc >> 1) ^ POLYNOMIAL;
                      else
                       crc >>= 1;
                     }
                  crc32_table[i] = crc;
                }
            }

     // Function to compute CRC32 without zlib
             uint32_t digest_c(const char *message) {
          init_crc32_table();  // Ensure lookup table is initialized

                  uint32_t crc = 0xFFFFFFFF; // Initial value

                                                                                                            while (*message) {
                                                                                                                         crc = (crc >> 8) ^ crc32_table[(crc ^ (uint8_t)(*message)) & 0xFF];
                                                                                                                                message++;
                                                                                                                                     }

                                                                                                                                         return crc ^ 0xFFFFFFFF; // Final XOR value
                                                                                                                                         }
