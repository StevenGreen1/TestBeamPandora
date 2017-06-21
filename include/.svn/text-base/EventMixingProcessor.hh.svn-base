#ifndef EventMixingProcessor_h
#define EventMixingProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <iostream>
#include <cmath>

#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/LCCollection.h>
#include <IO/LCReader.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <UTIL/LCRelationNavigator.h>

// CALICE includes
#include "AhcMapper.hh"
#include "MappedContainer.hh"

// gear
#include "gear/LayerLayout.h"

// 
#include "UTIL/CellIDDecoder.h"
#include "UTIL/CellIDEncoder.h"


using namespace lcio ;
using namespace marlin ;

/** EventMixingProcessor processor for overlaying background to each bunch crossing of a bunch train.
 *
 *  A more detailed description ....
 *
 *  @author P. Speckmayer CERN
 *
 *  @param a parameter
 *
 */


typedef std::map<int,CalorimeterHitImpl*> SortedCells;
typedef std::vector<std::string> StringVector;

class EventMixingProcessor : public Processor {  
     
public:
 
    virtual Processor*  newProcessor() { return new EventMixingProcessor ; }
 
 
    EventMixingProcessor() ;
 
    /** Called at the begin of the job before anything is read.
     * Use to initialize the processor, e.g. book histograms.
     */
    virtual void init() ;
 
    /**
     *  @brief  returns the processor name
     * 
     */
    virtual const std::string & name() const { return Processor::name() ; }  

    /**
     *  @brief  process the run header
     * 
     *  @param  run, the run header
     */
    virtual void processRunHeader( LCRunHeader* run ) ;

    /**
     *  @brief  process the event
     * 
     *  @param  evt, the event to be processed
     */
    virtual void processEvent( LCEvent * evt );

 
    /**
     *  @brief  is called at the end 
     * 
     */
    virtual void end() ;

protected:

    unsigned int _nRun ;
    unsigned int _nEvt ;

    LCReader* mixEventfileReader;

    std::string _ecalColSplitCaloHitsRelation;
    std::string _hcalColSplitCaloHitsRelation;
    std::string _tcColSplitCaloHitsRelation;

    StringVector _ecalSourceColNames;
    StringVector _hcalSourceColNames;
    StringVector _tcSourceColNames;

    std::string _ecalDestColName;
    std::string _hcalDestColName;
    std::string _tcDestColName;

    std::string _inputFileName;

    float _hcalCellSizeX;
    float _hcalCellSizeY;

    float _ecalCellSizeX;
    float _ecalCellSizeY;

    float _tcCellSizeX;
    float _tcCellSizeY;

    int _hcalMoveNCellsX;
    int _hcalMoveNCellsY;

    int _ecalMoveNCellsX;
    int _ecalMoveNCellsY;

    int _tcMoveNCellsX;
    int _tcMoveNCellsY;

    int _skipNEvents;

    std::string _mappingProcessorName;

    SortedCells sortedSplitCells;

    const CALICE::AhcMapper* _mapper;

//    inline long long cellID2long(int id0, int id1) { return ((long long) id0 << 32) | id1; };


    void combineCaloHits (EVENT::LCEvent* pLCEvent, EVENT::LCEvent* pLCEvent_add, StringVector& sourceCollectionNames, std::string& destCollectionName, 
                          std::string& relationCollectionName, float cell_size_x, float cell_size_y, 
                          int moveNCellsX, int moveNCellxY, const gear::LayerLayout &layerLayout, bool getCellSizesFromGear);
    
    void splitAndAddCells (LCCollection* sourceCollection, LCRelationNavigator& relationNavigator, float cellSizeU, float cellSizeV, 
                           int moveNCellsU, int moveNCellsV, UTIL::CellIDEncoder<CalorimeterHitImpl>& cellIdEncoder, const gear::LayerLayout &layerLayout, bool getCellSizesFromGear);

} ;

#endif
 
