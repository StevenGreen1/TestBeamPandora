/**
 *  @file   CalicePandora/src/CalicePandoraProcessor.cc
 * 
 *  @brief  Implementation of the pandora pfa new processor class.
 * 
 *  $Log: $
 */

#include "marlin/Global.h"
#include "marlin/Exceptions.h"

#include "gear/BField.h"

#include "Api/PandoraApi.h"

#include "CalicePseudoLayerCalculator.h"

#include "FineGranularityContent.h"


//#include "ExternalClusteringAlgorithm.h"
#include "PathLengthCalculator.h"
#include "CalicePandoraProcessor.h"
#include "SimpleBFieldCalculator.h"

#include <cstdlib>

CalicePandoraProcessor pandoraPFANewProcessor;

pandora::Pandora *CalicePandoraProcessor::m_pPandora = NULL;
EVENT::LCEvent *CalicePandoraProcessor::m_pLcioEvent = NULL;

//------------------------------------------------------------------------------------------------------------------------------------------

CalicePandoraProcessor::CalicePandoraProcessor() :
    Processor("CalicePandoraProcessor"),
    m_pGeometryCreator(NULL),
    m_pCaloHitCreator(NULL),
    m_pTrackCreator(NULL),
    m_pMCParticleCreator(NULL),
    m_pPfoCreator(NULL),
    m_nRun(0),
    m_nEvent(0)
{
    _description = "Pandora reconstructs clusters and particle flow objects";
    this->ProcessSteeringFile();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalicePandoraProcessor::init()
{
    try
    {
        streamlog_out(MESSAGE) << "CalicePandoraProcessor - Init" << std::endl;
        this->FinaliseSteeringParameters();

        streamlog_out(MESSAGE) << "CalicePandoraProcessor - Init - new pandora" << std::endl; // AAA0
        m_pPandora = new pandora::Pandora();
        streamlog_out(MESSAGE) << "CalicePandoraProcessor - Init - geometry creator" << std::endl; // AAA0
        m_pGeometryCreator = new GeometryCreator(m_geometryCreatorSettings);
        streamlog_out(MESSAGE) << "CalicePandoraProcessor - Init - calo hit creator" << std::endl; // AAA0
        m_pCaloHitCreator = new CaloHitCreator(m_caloHitCreatorSettings);
        streamlog_out(MESSAGE) << "CalicePandoraProcessor - Init - track creator" << std::endl; // AAA0
        m_pTrackCreator = new TrackCreator(m_trackCreatorSettings);
        streamlog_out(MESSAGE) << "CalicePandoraProcessor - Init - particle creator" << std::endl; // AAA0
        m_pMCParticleCreator = new MCParticleCreator(m_mcParticleCreatorSettings);
        streamlog_out(MESSAGE) << "CalicePandoraProcessor - Init - pfo creator" << std::endl; // AAA0
        m_pPfoCreator = new PfoCreator(m_pfoCreatorSettings);

        streamlog_out(MESSAGE) << "CalicePandoraProcessor - Init - RegisterUserComponents" << std::endl; // AAA0
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->RegisterUserComponents());
        streamlog_out(MESSAGE) << "CalicePandoraProcessor - Init - Create Geometry" << std::endl; // AAA1
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pGeometryCreator->CreateGeometry());
        streamlog_out(MESSAGE) << "CalicePandoraProcessor - Init - Read Settings" << std::endl; // AAA2
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*m_pPandora, m_settings.m_pandoraSettingsXmlFile));
    }
    catch (pandora::StatusCodeException &statusCodeException)
    {
        streamlog_out(ERROR) << "Failed to initialize Calice pandora: " << statusCodeException.ToString() << std::endl;
        throw statusCodeException;
    }
    catch (gear::Exception &exception)
    {
        streamlog_out(ERROR) << "Failed to initialize Calice pandora: gear exception " << exception.what() << std::endl;
        throw exception;
    }
    catch (...)
    {
        streamlog_out(ERROR) << "Failed to initialize Calice pandora: unrecognized exception" << std::endl;
        throw;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalicePandoraProcessor::processRunHeader(LCRunHeader *pLCRunHeader)
{
    m_detectorName = pLCRunHeader->getDetectorName();
    streamlog_out(MESSAGE) << "Detector Name " << m_detectorName << ", Run " << ++m_nRun <<  std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalicePandoraProcessor::processEvent(LCEvent *pLCEvent)
{
    static int eventCounter = 0;
    m_pLcioEvent = pLCEvent;

    if (eventCounter < m_settings.m_nEventsToSkip)
    {
        ++eventCounter;
        throw marlin::SkipEventException(this);
    }

    try
    {
        streamlog_out(MESSAGE) << "CalicePandoraProcessor, Run " << m_nRun << ", Event " << ++m_nEvent << "   lcio event number: " << pLCEvent->getEventNumber() << std::endl;

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pMCParticleCreator->CreateMCParticles(pLCEvent));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pTrackCreator->CreateTracks(pLCEvent));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pMCParticleCreator->CreateTrackToMCParticleRelationships(pLCEvent));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pCaloHitCreator->CreateCaloHits(pLCEvent));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pMCParticleCreator->CreateCaloHitToMCParticleRelationships(pLCEvent));

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pPandora));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pPfoCreator->CreateParticleFlowObjects(pLCEvent));

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pPandora));
        this->Reset();
    }
    catch (pandora::StatusCodeException &statusCodeException)
    {
        streamlog_out(ERROR) << "Calice pandora failed to process event: " << statusCodeException.ToString() << std::endl;
        throw statusCodeException;
    }
    catch (gear::Exception &exception)
    {
        streamlog_out(ERROR) << "Calice pandora failed to process event: gear exception " << exception.what() << std::endl;
        throw exception;
    }
    catch (EVENT::Exception &exception)
    {
        streamlog_out(ERROR) << "Calice pandora failed to process event: lcio exception " << exception.what() << std::endl;
        throw exception;
    }
    catch (...)
    {
        streamlog_out(ERROR) << "Calice pandora failed to process event: unrecognized exception" << std::endl;
        throw;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalicePandoraProcessor::check(LCEvent */*pLCEvent*/)
{
    streamlog_out(MESSAGE) << "CalicePandoraProcessor - Check" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalicePandoraProcessor::end()
{
    delete m_pPandora;
    delete m_pGeometryCreator;
    delete m_pCaloHitCreator;
    delete m_pTrackCreator;
    delete m_pMCParticleCreator;
    delete m_pPfoCreator;

    streamlog_out(MESSAGE) << "CalicePandoraProcessor - End" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode CalicePandoraProcessor::RegisterUserComponents() const
{
    // Register content from external pandora libraries
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetShowerProfileCalculator(*m_pPandora,
        new FineGranularityShowerProfileCalculator()));

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, FineGranularityContent::RegisterAlgorithms(*m_pPandora));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, FineGranularityContent::RegisterHelperFunctions(*m_pPandora));

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerCalculator(*m_pPandora, new CalicePseudoLayerCalculator()));

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetBFieldCalculator(*m_pPandora,
        new SimpleBFieldCalculator()));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void CalicePandoraProcessor::ProcessSteeringFile()
{
    registerProcessorParameter("PandoraSettingsXmlFile",
                            "The pandora settings xml file",
                            m_settings.m_pandoraSettingsXmlFile,
                            std::string());

    // Input collections
    registerProcessorParameter("TrackCollectionsX", 
                               "Names of the Track collection in X",
                               m_trackCreatorSettings.m_trackCollectionsX,
                               StringVector()); //std::string("TBTrackFEX"));

    registerProcessorParameter("TrackCollectionsY", 
                               "Names of the Track collection in Y",
                               m_trackCreatorSettings.m_trackCollectionsY,
                               StringVector()); // std::string("TBTrackFEY"));


//     registerInputCollections(LCIO::TRACK,
//                             "TrackCollections", 
//                             "Names of the Track collections used for clustering",
//                             m_trackCreatorSettings.m_trackCollections,
//                             StringVector());





    registerInputCollections(LCIO::CALORIMETERHIT,
                            "ECalCaloHitCollections", 
                            "Name of the ECAL calo hit collections",
                            m_caloHitCreatorSettings.m_eCalCaloHitCollections,
                            StringVector());

    registerInputCollections(LCIO::CALORIMETERHIT,
                            "HCalCaloHitCollections", 
                            "Name of the HCAL calo hit collections",
                            m_caloHitCreatorSettings.m_hCalCaloHitCollections,
                            StringVector());

    registerInputCollections(LCIO::CALORIMETERHIT,
                            "TcCaloHitCollections", 
                            "Name of the tail catcher calo hit collections",
                            m_caloHitCreatorSettings.m_tcCaloHitCollections,
                            StringVector());

    registerInputCollections(LCIO::MCPARTICLE,
                            "MCParticleCollections", 
                            "Name of mc particle collections",
                            m_mcParticleCreatorSettings.m_mcParticleCollections,
                            StringVector());

    registerInputCollections(LCIO::LCRELATION, 
                            "RelCaloHitCollections",
                            "SimCaloHit to CaloHit Relations Collection Name",
                            m_mcParticleCreatorSettings.m_lcCaloHitRelationCollections,
                            StringVector());

    registerInputCollections(LCIO::LCRELATION, 
                            "RelTrackCollections",
                            "Track to MCParticle Relations Collection Name",
                            m_mcParticleCreatorSettings.m_lcTrackRelationCollections,
                            StringVector());

    // Absorber properties ECal
    registerProcessorParameter("AbsorberRadiationLengthECal",
                            "The absorber radation length of the ECal",
                            m_geometryCreatorSettings.m_absorberRadiationLengthECal,
                            float(1.));

    registerProcessorParameter("AbsorberInteractionLengthECal",
                            "The absorber interaction length of the ECal",
                            m_geometryCreatorSettings.m_absorberInteractionLengthECal,
                            float(1.));

    // Absorber properties HCal
    registerProcessorParameter("AbsorberRadiationLengthHCal",
                            "The absorber radation length of the HCal",
                            m_geometryCreatorSettings.m_absorberRadiationLengthHCal,
                            float(1.));

    registerProcessorParameter("AbsorberInteractionLengthHCal",
                            "The absorber interaction length of the HCal",
                            m_geometryCreatorSettings.m_absorberInteractionLengthHCal,
                            float(1.));

    // Absorber properties tail catcher
    registerProcessorParameter("AbsorberRadiationLengthTc",
                            "The absorber radation length of the Tc",
                            m_geometryCreatorSettings.m_absorberRadiationLengthTc,
                            float(1.));

    registerProcessorParameter("AbsorberInteractionLengthTc",
                            "The absorber interaction length of the Tc",
                            m_geometryCreatorSettings.m_absorberInteractionLengthTc,
                            float(1.));

    // Name of PFO collection written by CalicePandora
    registerOutputCollection( LCIO::CLUSTER,
                              "ClusterCollectionName" , 
                              "Cluster Collection Name "  ,
                              m_pfoCreatorSettings.m_clusterCollectionName,
                              std::string("CalicePandoraClusters"));

    registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                              "PFOCollectionName" , 
                              "PFO Collection Name "  ,
                              m_pfoCreatorSettings.m_pfoCollectionName,
                              std::string("CalicePandoraPFOs"));

    // Calibration constants
    registerProcessorParameter("ECalToMipCalibration",
                            "The calibration from deposited ECal energy to mip",
                            m_caloHitCreatorSettings.m_eCalToMip,
                            float(1.));

    registerProcessorParameter("HCalToMipCalibration",
                            "The calibration from deposited HCal energy to mip",
                            m_caloHitCreatorSettings.m_hCalToMip,
                            float(1.));

    registerProcessorParameter("TcToMipCalibration",
                            "The calibration from deposited tail catcher energy to mip",
                            m_caloHitCreatorSettings.m_tcToMip,
                            float(1.));

    registerProcessorParameter("ECalMipThreshold",
                            "Threshold for creating calo hits in the ECal, units mip",
                            m_caloHitCreatorSettings.m_eCalMipThreshold,
                            float(0.));

    registerProcessorParameter("HCalMipThreshold",
                            "Threshold for creating calo hits in the HCal, units mip",
                            m_caloHitCreatorSettings.m_hCalMipThreshold,
                            float(0.));

    registerProcessorParameter("TcMipThreshold",
                            "Threshold for creating calo hits in the tail catcher, units mip",
                            m_caloHitCreatorSettings.m_tcMipThreshold,
                            float(0.));

    registerProcessorParameter("ECalToEMGeVCalibration",
                            "The calibration from deposited ECal energy to EM energy",
                            m_caloHitCreatorSettings.m_eCalToEMGeV,
                            float(1.));

    registerProcessorParameter("HCalToEMGeVCalibration",
                            "The calibration from deposited HCal energy to EM energy",
                            m_caloHitCreatorSettings.m_hCalToEMGeV,
                            float(1.));

    registerProcessorParameter("TcToEMGeVCalibration",
                            "The calibration from deposited tail catcher energy to EM energy",
                            m_caloHitCreatorSettings.m_tcToEMGeV,
                            float(1.));

    registerProcessorParameter("ECalToHadGeVCalibration",
                            "The calibration from deposited ECal energy to hadronic energy",
                            m_caloHitCreatorSettings.m_eCalToHadGeV,
                            float(1.));

    registerProcessorParameter("HCalToHadGeVCalibration",
                            "The calibration from deposited HCal energy to hadronic energy",
                            m_caloHitCreatorSettings.m_hCalToHadGeV,
                            float(1.));

    registerProcessorParameter("TcToHadGeVCalibration",
                            "The calibration from deposited tail catcher energy to hadronic energy",
                            m_caloHitCreatorSettings.m_tcToHadGeV,
                            float(1.));

    registerProcessorParameter("MaxHCalHitHadronicEnergy",
                            "The maximum hadronic energy allowed for a single hcal hit",
                            m_caloHitCreatorSettings.m_maxHCalHitHadronicEnergy,
                            float(1.));

    registerProcessorParameter("MaxTcHitHadronicEnergy",
                            "The maximum hadronic energy allowed for a single tail catcher hit",
                            m_caloHitCreatorSettings.m_maxTcHitHadronicEnergy,
                            float(1.));

    registerProcessorParameter("NOuterSamplingLayers",
                            "Number of layers from edge for hit to be flagged as an outer layer hit",
                            m_caloHitCreatorSettings.m_nOuterSamplingLayers,
                            int(3));

    registerProcessorParameter("LayersFromEdgeMaxRearDistance",
                            "Maximum number of layers from candidate outer layer hit to rear of detector",
                            m_caloHitCreatorSettings.m_layersFromEdgeMaxRearDistance,
                            float(250.f));

    // Average radiation length parameters
    registerProcessorParameter("AverageRadiationLengthECal",
                            "Average number of radiation lengths per mm in the ECal",
                            PathLengthCalculator::Settings::m_avgRadLengthECal,
                            float(0.13355)); // ILD_00, CLIC_ILD_01, CLIC_CDR value

    registerProcessorParameter("AverageRadiationLengthHCal",
                            "Average number of radiation lengths per mm in the HCal",
                            PathLengthCalculator::Settings::m_avgRadLengthHCal,
                            float(0.17392)); // CLIC_ILD_01, CLIC_CDR value. ILD_00 value is 0.043549

    // Average interaction length parameters
    registerProcessorParameter("AverageInteractionLengthECal",
                            "Average number of interaction lengths per mm in the ECal",
                            PathLengthCalculator::Settings::m_avgIntLengthECal,
                            float(0.0053555)); // ILD_00, CLIC_ILD_01, CLIC_CDR value

    registerProcessorParameter("AverageInteractionLengthHCal",
                            "Average number of interaction lengths per mm in the HCal",
                            PathLengthCalculator::Settings::m_avgIntLengthHCal,
                            float(0.0066047)); // CLIC_ILD_01, CLIC_CDR value. ILD_00 value is 0.0048187

    // track specifications
   registerProcessorParameter("TrackStartInZ",
                            "Z plane of virtual start of the track",
                            m_trackCreatorSettings.m_trackStartInZ,
                            float(-10000.f));

   registerProcessorParameter("TrackEndInZ",
                            "Z plane of virtual end of the track",
                            m_trackCreatorSettings.m_trackEndInZ,
                            float(0.f));

   registerProcessorParameter("TrackMomentumMagnitude",
                            "Momentum of the track (beam momentum)",
                            m_trackCreatorSettings.m_momentumMagnitude,
                            float(0.f));


    // Additional geometry parameters
    registerProcessorParameter("ECalSymmetryOrder",
                            "ECal symmetry order (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_eCalSymmetryOrder,
                            int(4));

    registerProcessorParameter("ECalPhiCoordinate",
                            "ECal phi coordinate (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_eCalPhiCoordinate,
                            float(0.));

    registerProcessorParameter("HCalSymmetryOrder",
                            "HCal symmetry order (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_hCalSymmetryOrder,
                            int(16));

    registerProcessorParameter("HCalPhiCoordinate",
                            "HCal phi coordinate (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_hCalPhiCoordinate,
                            float(0.));


    registerProcessorParameter("MarkTcHitsAsMuonHits",
                            "If true, mark tail catcher hits as pandora::MUON hits. If false, mark tail catcher hits as pandora::HCAL",
                            m_caloHitCreatorSettings.m_makeTcHitsMuonHits,
                            bool(false));


    // Number of events to skip
    registerProcessorParameter("NEventsToSkip",
                            "Number of events to skip at start of reconstruction job",
                            m_settings.m_nEventsToSkip,
                            int(0));


    // CALICE processors
    registerProcessorParameter( "MappingProcessorName" ,
                            "name of Ahcal MappingProcessor which takes care of the mapping",
                            m_caloHitCreatorSettings.m_mappingProcessorName,
                            std::string("AhcalMappingProcessor") ) ;

}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalicePandoraProcessor::FinaliseSteeringParameters()
{
    // ATTN: This function seems to be necessary for operations that cannot easily be performed at construction of the processor,
    // when the steering file is parsed e.g. the call to GEAR to get the inner bfield
    m_caloHitCreatorSettings.m_absorberRadiationLengthECal = m_geometryCreatorSettings.m_absorberRadiationLengthECal;
    m_caloHitCreatorSettings.m_absorberInteractionLengthECal = m_geometryCreatorSettings.m_absorberInteractionLengthECal;
    m_caloHitCreatorSettings.m_absorberRadiationLengthHCal = m_geometryCreatorSettings.m_absorberRadiationLengthHCal;
    m_caloHitCreatorSettings.m_absorberInteractionLengthHCal = m_geometryCreatorSettings.m_absorberInteractionLengthHCal;
    m_caloHitCreatorSettings.m_absorberRadiationLengthTc = m_geometryCreatorSettings.m_absorberRadiationLengthTc;
    m_caloHitCreatorSettings.m_absorberInteractionLengthTc = m_geometryCreatorSettings.m_absorberInteractionLengthTc;
//     m_caloHitCreatorSettings.m_hCalSymmetryOrder = m_geometryCreatorSettings.m_hCalSymmetryOrder;
//     m_caloHitCreatorSettings.m_hCalPhiCoordinate = m_geometryCreatorSettings.m_hCalPhiCoordinate;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalicePandoraProcessor::Reset()
{
    m_pCaloHitCreator->Reset();
    m_pTrackCreator->Reset();
}
