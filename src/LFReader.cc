#include <iostream>

#include "LFReader.h"
#include "LFFrames.h"

void printBits(size_t const size, void const * const ptr)
{
    unsigned char *b = (unsigned char*) ptr;
    unsigned char byte;
    int i, j;
    
    for (i = size-1; i >= 0; i--) {
        for (j = 7; j >= 0; j--) {
            byte = (b[i] >> j) & 1;
            printf("%u", byte);
        }
    }
    puts("");
}

LFReader::LFReader( std::map<int,int> boards, bool isPileUp )
{

  fBoards = boards;

  // To initialize the flags
  for( auto it : fBoards ){
    fROOTInitialized[it.first] = false;
  }
    
}

bool LFReader::InitializeROOT( uint32_t& bIdx, uint32_t& aggSize ){

  // Initializing the histograms
  //std::cout << "Initializing ROOT objects for board " << bIdx << std::endl;

  // Creating ROOT file
  fFileOut = new TFile( outFile.c_str( ), "RECREATE" ); 

  if( fBoards.find(bIdx) == fBoards.end() ) {
    std::cout << "Board " << bIdx << " not found." << std::endl;
    return false;
  }

  for(UInt_t j = 0 ; j < fBoards[bIdx]; j++){
    if( aggSize == sizeof(caenEventPHA_t) ){
      fEnergy[bIdx].push_back(new TH1F(Form("hEnergy_%i_%i",bIdx,j),
			        Form("Energy Board %i Channel %i",bIdx,j),
			        32768,0,32768));
      fEnergy[bIdx].back()->SetXTitle("Energy [channels]");
      fEnergy[bIdx].back()->SetYTitle("Counts / channel");
    }
    else if( aggSize == sizeof(caenEventPSD_t) ){
      fQshort[bIdx].push_back(new TH1F(Form("hQshort_%i_%i",bIdx,j),
			        Form("Qshort Board %i Channel %i",bIdx,j),
			        32768,0,32768));
      fQlong[bIdx].push_back(new TH1F(Form("hQlong_%i_%i",bIdx,j),
			        Form("Qlong Board %i Channel %i",bIdx,j),
			        65535,0,65535));
      fQshort[bIdx].back()->SetXTitle("Qshort [channels]");
      fQshort[bIdx].back()->SetYTitle("Counts / channel");
      fQlong[bIdx].back()->SetXTitle("Qlong [channels]");
      fQlong[bIdx].back()->SetYTitle("Counts / channel");
    }

  }

  if( aggSize == sizeof(caenEventPHA_t) ){
    fTree[bIdx] = new TTree(Form("Board %i",bIdx),Form("Board %i",bIdx));
    fTree[bIdx]->Branch( "pu"       ,        &fPu[bIdx],        "pu/O" );
    fTree[bIdx]->Branch( "satu"     ,      &fSatu[bIdx],      "satu/O" );
    fTree[bIdx]->Branch( "lost"     ,      &fLost[bIdx],      "lost/O" );
    fTree[bIdx]->Branch( "channel"  ,   &fChannel[bIdx],   "channel/s" );
    fTree[bIdx]->Branch( "cfd"      ,       &fCFD[bIdx],       "cfd/s" );
    fTree[bIdx]->Branch( "timeStamp", &fTimeStamp[bIdx], "timeStamp/l" );
    fTree[bIdx]->Branch( "energy"   ,  &fEnergies[bIdx],    "energy/s" );
    fTree[bIdx]->SetMaxVirtualSize(1000000000LL);
  }
  else if( aggSize == sizeof(caenEventPSD_t) ){
    fTree[bIdx] = new TTree(Form("Board %i",bIdx),Form("Board %i",bIdx));
    fTree[bIdx]->Branch( "pu"       ,        &fPu[bIdx],        "pu/O" );
    fTree[bIdx]->Branch( "satu"     ,      &fSatu[bIdx],      "satu/O" );
    fTree[bIdx]->Branch( "lost"     ,      &fLost[bIdx],      "lost/O" );
    fTree[bIdx]->Branch( "channel"  ,   &fChannel[bIdx],   "channel/s" );
    fTree[bIdx]->Branch( "cfd"      ,       &fCFD[bIdx],       "cfd/s" );
    fTree[bIdx]->Branch( "timeStamp", &fTimeStamp[bIdx], "timeStamp/l" );
    fTree[bIdx]->Branch( "qShort"   ,    &fQShort[bIdx],    "qShort/s" );
    fTree[bIdx]->Branch( "qLong"    ,     &fQLong[bIdx],     "qLong/s" );
    fTree[bIdx]->SetMaxVirtualSize(1000000000LL);
  }

  fROOTInitialized[bIdx] = true;

  return true;

}

void LFReader::UnpackPHA( char* buffer, uint32_t& offset, uint32_t& boardId, uint16_t& chan ){

  memcpy( &fDataPHA, buffer+offset+sizeof(dataKey), sizeof(fDataPHA) );

  fEnergies[boardId] = fDataPHA.energy;
  fEnergy[boardId][chan]->Fill( fDataPHA.energy );

  fCFD[boardId] = fDataPHA.cfd;

}

void LFReader::UnpackPSD( char* buffer, uint32_t& offset, uint32_t& boardId, uint16_t& chan ){

  memcpy( &fDataPSD, buffer+offset+sizeof(dataKey), sizeof(fDataPSD) );

  fQLong[boardId] = fDataPSD.qlong;
  fQlong[boardId][chan]->Fill( fDataPSD.qlong );

  fQShort[boardId] = fDataPSD.qshort;
  fQshort[boardId][chan]->Fill( fDataPSD.qshort );

  fCFD[boardId] = fDataPSD.cfd;

}

bool LFReader::Read( std::string in, std::string out ){

  outFile = out;

  std::ifstream input( in, std::ios::binary );
  if( !input.is_open( ) ){    
    std::cerr << "Input File not open please check the input file name!" << std::endl;
    return 0;
  }

  uint64_t length;
  input.seekg(0, std::ios::end);
  length = (uint64_t)input.tellg();
  
  uint32_t aggSize, boardId;

  uint64_t bufferSize = 1e8;
  char* buffer = new char[bufferSize];

  uint32_t ts       = 0;
  uint32_t tsoffset = 0;
  
  double progress;
  uint64_t pos = 0;
  ULong_t i, posBar, barWidth = 70;
  while( pos < length ){

    input.seekg( pos, std::ios::beg );
    input.read( buffer, bufferSize );

    uint32_t offset = 0;
    while( offset < bufferSize ){

      /* Progress Bar */ 
      progress = (double)(offset+pos)/(double)length;
      if( (int)(progress*100) % 1 == 0 ){
        posBar = barWidth * progress;
        std::cout << "[";
        for( int k = 0; k < barWidth; ++k ){
	        if( k < posBar ) std::cout << "=";
	        else if( k == posBar ) std::cout << ">";
	        else std::cout << " ";
        }
        std::cout << "] " << int( progress * 100 ) << "%\r";
        std::cout.flush( );
      }

      if( offset + pos >= length ){
        pos += offset;
        break;
      }

      memcpy( &fKey, buffer+offset, sizeof( fKey ) );
      
      aggSize = fKey.GetBytes( );
      boardId = fKey.GetBoard( );

      if( aggSize + offset >= bufferSize ){
        pos += offset;
        break;
      }
      
      fChannel[boardId]   = fKey.GetChannel( );
      fTimeStamp[boardId] = fKey.GetTstamp( );
      fPu[boardId]        = fKey.GetPu( );
      fSatu[boardId]      = fKey.GetSatu( );
      fLost[boardId]      = fKey.GetLost( );
      
      if( !fROOTInitialized[boardId] ){
        if( !InitializeROOT( boardId, aggSize ) ){
          return false;
        }
      }

      if( fKey.IsIdle() ){
	      offset += aggSize;
	      continue;
      }
      
      if( aggSize == sizeof(caenEventPHA_t) )
	      UnpackPHA( buffer, offset, boardId, fChannel[boardId] );
      else if( aggSize == sizeof(caenEventPSD_t) )
	      UnpackPSD( buffer, offset, boardId, fChannel[boardId] );
      
      offset += aggSize;
      fTree[boardId]->Fill( );

    }

  }

  std::cout << std::endl;

  return true;

}

void LFReader::Write( ){

  std::cout << "Saving data to ROOT file..." << std::endl;
  
  fFileOut->cd( );

  if( !fEnergy.empty( ) ){
    fFileOut->mkdir( "Energies" );
    for( auto it1 : fEnergy ){
      fFileOut->cd( "Energies" );
      for( auto it2 : it1.second ){
        it2->Write( );
      }
    }
  }

  if( !fQlong.empty( ) ){
    fFileOut->mkdir( "Qlong" );
    for( auto it1 : fQlong ){
      fFileOut->cd( "Qlong" );
      for( auto it2 : it1.second ){
        it2->Write( );
      }
    }
  }

  if( !fQshort.empty( ) ){
    fFileOut->mkdir( "Qshort" );
    for( auto it1 : fQshort ){
      fFileOut->cd( "Qshort" );
      for( auto it2 : it1.second ){
        it2->Write( );
      }
    }
  }

  fFileOut->cd( );
  for( auto it : fTree ){
    it.second->Write( 0, TObject::kOverwrite );
  }

  fFileOut->Close( );

  std::cout << "Saved data to ROOT file." << std::endl << std::endl;
  

}

