#ifndef __OFDM__Transmitter__
#define __OFDM__Transmitter__

#include <fstream>
#include "ofdm_params.h"

class Transmitter {
public:
    Transmitter();
    
    bool ReadDataFromFile( const char* filename );
    bool WriteDataToFile( const char* filename );
    void GenerateHeader( );
    void Modulate( char* pData, unsigned int uDataLen, unsigned int uIfftSize, unsigned int uCarriers, unsigned int uConjCarriers, unsigned int uCarrierCount, unsigned int uSymbSize, unsigned int uGuardTime );
    void TransmitData( );
    
private:
    unsigned char* m_pData;
    unsigned long m_lDataLength;
    
};

#endif /* defined(__OFDM__Transmitter__) */
