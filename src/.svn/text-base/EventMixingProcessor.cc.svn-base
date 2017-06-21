#include "marlin/Global.h"

#include "EventMixingProcessor.hh"
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>

#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCFlagImpl.h>

#include "MappingProcessor.hh"
#include <DecoderSet.hh>

#include "gear/GearParameters.h"
#include "gear/CalorimeterParameters.h"
#include "gear/GEAR.h"
#include "gear/LayerLayout.h"


// ----- include for verbosity dependend logging ---------

using namespace lcio ;
using namespace marlin ;


EventMixingProcessor aEventMixingProcessor ;

EventMixingProcessor::EventMixingProcessor() : Processor("EventMixingProcessor") {
 
  // modify processor description
  _description = "Processeor to mix in Calice events from another file; a shift in X and Y can be applied" ;

 StringVec files ;


 // CALICE processors
 registerProcessorParameter( "MappingProcessorName" ,
                             "Name of Ahcal MappingProcessor which takes care of the mapping",
                             _mappingProcessorName,
                             std::string("AhcalMappingProcessor") ) ;
 // split cell relations
  registerProcessorParameter("ECalSplitRelations",
                             "Name of the collection which stores the relations from the new (split and potentially moved) ECal cells to the original cells"  ,
                             _ecalColSplitCaloHitsRelation,
                             std::string("ECalSplitRelations")) ;

  registerProcessorParameter("HCalSplitRelations",
                             "Name of the collection which stores the relations from the new (split and potentially moved) HCal cells to the original cells"  ,
                             _hcalColSplitCaloHitsRelation,
                             std::string("HCalSplitRelations")) ;

  registerProcessorParameter("TCSplitRelations",
                             "Name of the collection which stores the relations from the new (split and potentially moved) tail catcher cells to the original cells"  ,
                             _tcColSplitCaloHitsRelation,
                             std::string("TCSplitRelations")) ;



  // source collections
  registerProcessorParameter("ECalSourceCollections",
                             "Name of the source collections for the ECal"  ,
                             _ecalSourceColNames,
                             StringVector()) ;

  registerProcessorParameter("HCalSourceCollections",
                             "Name of the source collections for the HCal"  ,
                             _hcalSourceColNames,
                             StringVector()) ;

  registerProcessorParameter("TCSourceCollections",
                             "Name of the source collections for the tail catcher"  ,
                             _tcSourceColNames,
                             StringVector()) ;


  // destination collections
  registerProcessorParameter("ECalDestinationCollection",
                             "Name of the source collection for the ECal"  ,
                             _ecalDestColName,
                             std::string("SplitECalCaloHits")) ;

  registerProcessorParameter("HCalDestinationCollection",
                             "Name of the source collection for the HCal"  ,
                             _hcalDestColName,
                             std::string("SplitHCalCaloHits")) ;

  registerProcessorParameter("TCDestinationCollection",
                             "Name of the source collection for the tail catcher"  ,
                             _tcDestColName,
                             std::string("SplitTCCaloHits")) ;



  // mix file name
  registerProcessorParameter( "MixFileName" ,
                              "Name of the lcio input file with the events to mix."  ,
                              _inputFileName,
                              std::string("")) ;



  // split cell sizes HCal
  registerProcessorParameter("HCalCellSizeX",
                             "New cell size X of HCal cells"  ,
                             _hcalCellSizeX,
                             float(3.f)) ;

  registerProcessorParameter("HCalCellSizeY",
                             "New cell size Y of HCal cells"  ,
                             _hcalCellSizeY,
                             float(3.f)) ;

  // split cell sizes ECal
  registerProcessorParameter("ECalCellSizeX",
                             "New cell size X of ECal cells"  ,
                             _ecalCellSizeX,
                             float(1.f)) ;

  registerProcessorParameter("ECalCellSizeY",
                             "New cell size Y of ECal cells"  ,
                             _ecalCellSizeY,
                             float(1.f)) ;

  // split cell sizes tail catcher
  registerProcessorParameter("TCCellSizeX",
                             "New cell size X of tail catcher cells"  ,
                             _tcCellSizeX,
                             float(0.f)) ; // 0. --> don't split cells

  registerProcessorParameter("TCCellSizeY",
                             "New cell size Y of tail catcher cells"  ,
                             _tcCellSizeY,
                             float(0.f)) ; // 0. --> don't split cells





  // move cells HCal
  //
  // BE CAREFUL when moving cells by negative values. The outermost cells will hit easily the lower limit of the BitField of which the value cannot be lower than 1 !!!!!!!!!
  //
  registerProcessorParameter("MoveHCalNCellsX",
                             "Move HCal hits by N cells in the X direction"  ,
                             _hcalMoveNCellsX,
                             int(0.f)) ;

  registerProcessorParameter("MoveHCalNCellsY",
                             "Move HCal hits by N cells in the X direction"  ,
                             _hcalMoveNCellsY,
                             int(0.f)) ;

  // move cells ECal
  registerProcessorParameter("MoveECalNCellsX",
                             "Move ECal hits by N cells in the X direction"  ,
                             _ecalMoveNCellsX,
                             int(0.f)) ;

  registerProcessorParameter("MoveECalNCellsY",
                             "Move ECal hits by N cells in the Y direction"  ,
                             _ecalMoveNCellsY,
                             int(0.f)) ;

  // move cells tail catcher
  registerProcessorParameter("MoveTCNCellsX",
                             "Move tail catcher hits by N cells in the X direction"  ,
                             _tcMoveNCellsX,
                             int(0.f)) ;

  registerProcessorParameter("MoveTCNCellsY",
                             "Move tail catcher hits by N cells in the Y direction"  ,
                             _tcMoveNCellsY,
                             int(0.f)) ;


  registerProcessorParameter("SkipNEvents",
                             "Skip given number of events"  ,
                             _skipNEvents,
                             int(0.f)) ;


}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void EventMixingProcessor::init() {
    streamlog_out(DEBUG) << " init called  " << std::endl ;
    printParameters() ;
 
    mixEventfileReader = LCFactory::getInstance()->createLCReader() ;
    mixEventfileReader->open (_inputFileName);
    mixEventfileReader->skipNEvents(_skipNEvents);
 
    _nRun = 0 ;
    _nEvt = 0 ;

    _mapper = dynamic_cast<const CALICE::AhcMapper*> ( CALICE::MappingProcessor::getMapper(_mappingProcessorName) );
    if (!_mapper) {
        streamlog_out(ERROR) << "Cannot obtain AhcMapper from MappingProcessor " << _mappingProcessorName
                             <<". Mapper not present or wrong type." << std::endl;
    }
}


void EventMixingProcessor::processRunHeader( LCRunHeader* run) 
{
    _nRun++ ;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void EventMixingProcessor::processEvent( LCEvent * evt )
{
  LCEvent * mixEvt = mixEventfileReader->readNextEvent( LCIO::UPDATE ) ;

  if (mixEvt==0)
  {
      mixEventfileReader->close();
      return;
  }


  if (!_ecalDestColName.empty() && !_ecalColSplitCaloHitsRelation.empty())
  {
      const gear::LayerLayout &ecalLayerLayout(marlin::Global::GEAR->getEcalEndcapParameters().getLayerLayout()); 
      combineCaloHits (evt, mixEvt, _ecalSourceColNames, _ecalDestColName, _ecalColSplitCaloHitsRelation, _ecalCellSizeX, _ecalCellSizeY,
                       _ecalMoveNCellsX, _ecalMoveNCellsY, ecalLayerLayout, true);
  }
                   
  if (!_hcalDestColName.empty() && !_hcalColSplitCaloHitsRelation.empty())
  {
      const gear::LayerLayout &hcalLayerLayout(marlin::Global::GEAR->getHcalEndcapParameters().getLayerLayout()); 
      combineCaloHits (evt, mixEvt, _hcalSourceColNames, _hcalDestColName, _hcalColSplitCaloHitsRelation, _hcalCellSizeX, _hcalCellSizeY,
                       _hcalMoveNCellsX, _hcalMoveNCellsY, hcalLayerLayout, false);
  }

  if (!_tcDestColName.empty() && !_tcColSplitCaloHitsRelation.empty())
  {
      const gear::LayerLayout &tcLayerLayout(marlin::Global::GEAR->getYokeEndcapParameters().getLayerLayout()); 
      combineCaloHits (evt, mixEvt, _tcSourceColNames, _tcDestColName, _tcColSplitCaloHitsRelation, _tcCellSizeX, _tcCellSizeY,
                       _tcMoveNCellsX, _tcMoveNCellsY, tcLayerLayout, true);
  }                   

  ++_nEvt;  
}







/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void EventMixingProcessor::combineCaloHits (EVENT::LCEvent* pLCEvent, EVENT::LCEvent* pLCEvent_add, StringVector& sourceCollectionNames, std::string& destCollectionName, 
                                            std::string& relationCollectionName,
                                            float cell_size_x, float cell_size_y, int moveNCellsX, int moveNCellsY, const gear::LayerLayout &layerLayout, bool getCellSizesFromGear)
{
//    std::cout << "Event  : " << pLCEvent->getEventNumber() << std::endl;

    LCRelationNavigator relationNavigator(LCIO::CALORIMETERHIT, LCIO::CALORIMETERHIT);


    std::string cellEncoding("");
    // get the cell encoding
    for (StringVector::iterator itCollName=sourceCollectionNames.begin(), itCollNameEnd=sourceCollectionNames.end(); itCollName!=itCollNameEnd; ++itCollName)
    {
        std::string sourceName=(*itCollName);
        try
        {
            EVENT::LCCollection *pCaloHitCollection = pLCEvent->getCollection(sourceName);

            // cell encoding of source collection
            if (cellEncoding.empty())
            {
                UTIL::CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);
                cellEncoding = pCaloHitCollection->getParameters().getStringVal("CellIDEncoding");

                _mapper->getDecoder()->setCellIDEncoding(pCaloHitCollection->getParameters().getStringVal(LCIO::CellIDEncoding));

                break;
            }
        }
        catch (EVENT::Exception &exception)
        {
            std::cout << "Failed to extract collection. " << exception.what() << std::endl;
        }
    }


    // could not extract cell encoding, collection is probably not existing in this event
    if (cellEncoding.empty())
        return;



    // grab or create the destination collection
    EVENT::LCCollection *pCaloHitCollectionDestination = NULL;
    try
    {
        pCaloHitCollectionDestination = pLCEvent->getCollection(destCollectionName);
    }
    catch (EVENT::Exception &exception)
    {
        streamlog_out(DEBUG) << "Failed to extract destination calo hit collection --> create it.    " << exception.what() << std::endl;
    }


    if (!pCaloHitCollectionDestination)
    {
        pCaloHitCollectionDestination = new LCCollectionVec(LCIO::CALORIMETERHIT);

        /*we want to save the position, and this can only be done if the CHBIT_LONG bit is set*/
        pCaloHitCollectionDestination->setFlag(pCaloHitCollectionDestination->getFlag() | 1 << LCIO::CHBIT_LONG );
    }


    // check if the destination collection is now really there
    if (!pCaloHitCollectionDestination)
        throw EVENT::Exception("Destination collection could not be created");

    // cell encoding of destination collection
    UTIL::CellIDEncoder<CalorimeterHitImpl> cellIdEncoder( cellEncoding, pCaloHitCollectionDestination ) ;

    // add cells from the pLCEvent
    for (StringVector::iterator itCollName=sourceCollectionNames.begin(), itCollNameEnd=sourceCollectionNames.end(); itCollName!=itCollNameEnd; ++itCollName)
    {
        std::string sourceName=(*itCollName);
        try
        {
            EVENT::LCCollection *pCaloHitCollection = pLCEvent->getCollection(sourceName);
            streamlog_out(DEBUG) << "retrieve base event" << std::endl;
            splitAndAddCells (pCaloHitCollection, relationNavigator, cell_size_x, cell_size_y, 0, 0, cellIdEncoder, layerLayout, getCellSizesFromGear);
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(MESSAGE) << "Failed to extract collection from base event. " << exception.what() << std::endl;
        }
    }


    // add cells from the pLCEvent_add
    for (StringVector::iterator itCollName=sourceCollectionNames.begin(), itCollNameEnd=sourceCollectionNames.end(); itCollName!=itCollNameEnd; ++itCollName)
    {
        std::string sourceName=(*itCollName);
        try
        {
            EVENT::LCCollection *pCaloHitCollection = pLCEvent_add->getCollection(sourceName);
            streamlog_out(DEBUG) << "retrieve additional event" << std::endl;
            splitAndAddCells (pCaloHitCollection, relationNavigator, cell_size_x, cell_size_y, moveNCellsX, moveNCellsY, cellIdEncoder, layerLayout, getCellSizesFromGear);
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(WARNING) << "Failed to extract collection from added event. " << exception.what() << std::endl;
        }
    }

    // write hits into destination collection
    for (SortedCells::iterator itCell=sortedSplitCells.begin(), itCellEnd=sortedSplitCells.end(); itCell != itCellEnd; ++itCell)
    {
        CalorimeterHitImpl* pCaloHit=itCell->second;
        pCaloHitCollectionDestination->addElement(pCaloHit);
    }
 
    EVENT::LCParameters &collParameters = pCaloHitCollectionDestination->parameters ();
    collParameters.setValue("UseGearForHCalCellSize", int(1));

    pLCEvent->addCollection(pCaloHitCollectionDestination, destCollectionName.c_str());
    pLCEvent->addCollection(relationNavigator.createLCCollection(), relationCollectionName);

    sortedSplitCells.clear();
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void EventMixingProcessor::splitAndAddCells (LCCollection* sourceCollection, LCRelationNavigator& relationNavigator, 
                                             float cellSizeU, float cellSizeV, int moveNCellsU, int moveNCellsV,  UTIL::CellIDEncoder<CalorimeterHitImpl>& cellIdEncoder, 
                                             const gear::LayerLayout &layerLayout, bool getCellSizesFromGear)
{
    const int nElements(sourceCollection->getNumberOfElements());

    if (0 == nElements)
        return;

    // cell encoding of source collection
    UTIL::CellIDDecoder<CalorimeterHit> cellIdDecoder(sourceCollection);
    std::string cellEncoding = sourceCollection->getParameters().getStringVal("CellIDEncoding");


    CalorimeterHitImpl* temporaryCaloHit = new CalorimeterHitImpl();


    //in loop
    for (int i = 0; i < nElements; ++i)
    {
        try
        {
            EVENT::CalorimeterHit *pCaloHit = dynamic_cast<CalorimeterHit*>(sourceCollection->getElementAt(i));

            const float *pCaloHitPosition(pCaloHit->getPosition());
            int cellID=pCaloHit->getCellID0();

            int cellIdI = cellIdDecoder(pCaloHit)["I"];
            int cellIdJ = cellIdDecoder(pCaloHit)["J"];
            int cellIdK = cellIdDecoder(pCaloHit)["K-1"]; // check if "K-1"+1


            float inputCellSizeI=0.f;
            float inputCellSizeJ=0.f;


            // get cell dimensions in U and V
            if (getCellSizesFromGear)
            {
                inputCellSizeI = layerLayout.getCellSize0 (cellIdK); // get the cell size I for layer K
                inputCellSizeJ = layerLayout.getCellSize1 (cellIdK); // get the cell size J for layer K
            }
            else
            {
                inputCellSizeI = _mapper->getISizeFromCellID (cellID)*10.0;
                inputCellSizeJ = _mapper->getJSizeFromCellID (cellID)*10.0;
            }


            // if cell sizes set to 0. then don't split the cells
            if (cellSizeU<=0.f)
                cellSizeU=inputCellSizeI;

            if (cellSizeV<=0.f)
                cellSizeV=inputCellSizeJ;

            // divide cell
            int divisionsU = int(floor(0.5 +inputCellSizeI/(cellSizeU*10.f))); 
            int divisionsV = int(floor(0.5 +inputCellSizeJ/(cellSizeV*10.f))); 

            int numberCells = divisionsU*divisionsV;

            float inputCellEnergy = pCaloHit->getEnergy();
            float splitEnergy     = inputCellEnergy/numberCells;

            // create the new cells
            for (int iCellU=0; iCellU<divisionsU; ++iCellU) 
            {
                for (int iCellV=0; iCellV<divisionsV; ++iCellV) 
                {

                    // calculate cell IDs
                    int cellIdSplitI = cellIdI +int(moveNCellsU*cellSizeU/1.0)+iCellU; // Calice cell ID unit is 1 cm
                    int cellIdSplitJ = cellIdJ +int(moveNCellsV*cellSizeV/1.0)+iCellV; // Calice cell ID unit is 1 cm
                    int cellIdSplitK = cellIdK; // the layer stays the same


                    cellIdEncoder["I"]   = cellIdSplitI ;
                    cellIdEncoder["J"]   = cellIdSplitJ ;
                    cellIdEncoder["K-1"] = cellIdSplitK ;


                    // set the cell Id of the new calohit
                    cellIdEncoder.setCellID (temporaryCaloHit);
                    int cellIdSplitCell=temporaryCaloHit->getCellID0();

//                    std::cout << "original cellID  "<< cellID << "  icellv" << iCellV << "  I " << cellIdI << "  J " << cellIdJ << "  K " << cellIdK << "   divU " << divisionsU << " divV " << divisionsV << "  cidsplitI " << cellIdSplitI << "  cellIdSplitJ " << cellIdSplitJ << "   cellIdsplitcell " << cellIdSplitCell << std::endl;


                    CalorimeterHitImpl* splitCaloHit = NULL;
                    // check if cell with same cell id is already present
                    SortedCells::iterator itCell = sortedSplitCells.find(cellIdSplitCell);
                    if (itCell!=sortedSplitCells.end())
                    {

                        // there is already a cell with this ID, add the splitEnergy to that cell
                        splitCaloHit = itCell->second;
                        float energy = splitCaloHit->getEnergy();
                        splitCaloHit->setEnergy(energy+splitEnergy);
                    }
                    else
                    {
                        // no cell with this ID present, create a new one
                        splitCaloHit = new CalorimeterHitImpl();

                        // calculate position of the new cell
                        float hitPosition[3];
                        hitPosition[0]=pCaloHitPosition[0]+float(moveNCellsU+iCellU)*cellSizeU*10.f;
                        hitPosition[1]=pCaloHitPosition[1]+float(moveNCellsV+iCellV)*cellSizeV*10.f;
                        hitPosition[2]=pCaloHitPosition[2];

                        splitCaloHit->setEnergy(splitEnergy);
                        splitCaloHit->setPosition(hitPosition);

                        // set the cell Id of the new calohit
                        cellIdEncoder.setCellID (splitCaloHit);
                        int cellId=splitCaloHit->getCellID0();

                        // add cell to sortedSplitCells map
                        sortedSplitCells.insert(std::make_pair(cellId,splitCaloHit));
                    }
                    // add a relation from the new cell to the old cell. Take as weight the splitenergy 
                    relationNavigator.addRelation(splitCaloHit, pCaloHit, splitEnergy);
        
                }
            }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(ERROR) << "Failed to extract calo hit: " << exception.what() << std::endl;
        }
    }

    delete temporaryCaloHit;
}






/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void EventMixingProcessor::end(){
    mixEventfileReader->close();
}

  
