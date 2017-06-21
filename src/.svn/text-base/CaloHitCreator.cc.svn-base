/**
 *  @file   MarlinPandora/src/CaloHitCreator.cc
 * 
 *  @brief  Implementation of the calo hit creator class.
 * 
 *  $Log: $
 */

#include "marlin/Global.h"
#include "marlin/Processor.h"

#include "gear/GearParameters.h"
#include "gear/CalorimeterParameters.h"
#include "gear/GEAR.h"
#include "gear/LayerLayout.h"

#include "UTIL/CellIDDecoder.h"

#include "MappingProcessor.hh"

#include "CaloHitCreator.h"
#include "PathLengthCalculator.h"
#include "CalicePandoraProcessor.h"

#include <cmath>
#include <limits>

CalorimeterHitVector CaloHitCreator::m_calorimeterHitVector;

CaloHitCreator::CaloHitCreator(const Settings &settings) :
    m_settings(settings),
    m_pPandora(CalicePandoraProcessor::GetPandora()),
    m_pPathLengthCalculator(PathLengthCalculator::GetInstance()),
    m_eCalZ(0),
    m_eCalR(0),
    m_eCalPhi0(0),
    m_eCalSymmetry(0),
    m_hCalZ(0),
    m_hCalR(0),
    m_hCalPhi0(0),
    m_hCalSymmetry(0),
    m_tcZ(0),
    m_tcR(0),
    m_tcPhi0(0),
    m_tcSymmetry(0),
    m_hasECal(true),
    m_hasHCal(true),
    m_hasTc(true),
    m_eCalLayerThickness(0),
    m_hCalLayerThickness(0),
    m_tcLayerThickness(0)
{
    try
    {
        m_eCalZ = marlin::Global::GEAR->getEcalEndcapParameters().getExtent()[3];
        m_eCalR = marlin::Global::GEAR->getEcalEndcapParameters().getExtent()[1];
        m_eCalPhi0 = marlin::Global::GEAR->getEcalEndcapParameters().getPhi0();
        m_eCalSymmetry = marlin::Global::GEAR->getEcalEndcapParameters().getSymmetryOrder();
        const gear::LayerLayout &eCalLayerLayout(marlin::Global::GEAR->getEcalEndcapParameters().getLayerLayout()); 
        m_eCalLayerThickness = eCalLayerLayout.getThickness(eCalLayerLayout.getNLayers() - 1);
        if (0.f == m_eCalLayerThickness)
        {
            streamlog_out(ERROR) << "CaloHitCreator/ m_eCalLayerThickness==0" << std::endl;
            throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
        }
    }
    catch (gear::UnknownParameterException &unknownParameterException)
    {
        streamlog_out(MESSAGE) << "CaloHitCreator/No ECal found in GEAR file : " << unknownParameterException.what() << std::endl;
        m_hasECal = false;
    }
    
    try
    {
        m_hCalZ = marlin::Global::GEAR->getHcalEndcapParameters().getExtent()[3];
        m_hCalR = marlin::Global::GEAR->getHcalEndcapParameters().getExtent()[1];
//         m_hCalPhi0 = marlin::Global::GEAR->getHcalEndcapParameters().getIntVal("Hcal_outer_polygon_phi0");
//         m_hCalSymmetry = marlin::Global::GEAR->getHcalEndcapParameters().getIntVal("Hcal_outer_polygon_order");
        m_hCalPhi0 = 0.f;
        m_hCalSymmetry = 4;
        const gear::LayerLayout &hCalLayerLayout(marlin::Global::GEAR->getHcalEndcapParameters().getLayerLayout()); 
        m_hCalLayerThickness = hCalLayerLayout.getThickness(hCalLayerLayout.getNLayers() - 1);
        if (0.f == m_hCalLayerThickness)
        {
            streamlog_out(ERROR) << "CaloHitCreator/ m_hCalLayerThickness==0" << std::endl;
            throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
        }
    }
    catch (gear::UnknownParameterException &unknownParameterException)
    {
        streamlog_out(MESSAGE) << "CaloHitCreator/No HCal found in GEAR file : " << unknownParameterException.what() << std::endl;
        m_hasHCal = false;
    }

    try
    {
        m_tcZ = marlin::Global::GEAR->getYokeEndcapParameters().getExtent()[3];
        m_tcR = marlin::Global::GEAR->getYokeEndcapParameters().getExtent()[1];
//         m_hCalPhi0 = marlin::Global::GEAR->getHcalEndcapParameters().getIntVal("Hcal_outer_polygon_phi0");
//         m_hCalSymmetry = marlin::Global::GEAR->getHcalEndcapParameters().getIntVal("Hcal_outer_polygon_order");
        m_tcPhi0 = 0.f;
        m_tcSymmetry = 4;
        const gear::LayerLayout &tcLayerLayout(marlin::Global::GEAR->getYokeEndcapParameters().getLayerLayout()); 
        m_tcLayerThickness = tcLayerLayout.getThickness(tcLayerLayout.getNLayers() - 1);
        if (0.f == m_tcLayerThickness)
        {
            streamlog_out(ERROR) << "CaloHitCreator/ m_tcLayerThickness==0" << std::endl;
            throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
        }
    }
    catch (gear::UnknownParameterException &unknownParameterException)
    {
        streamlog_out(MESSAGE) << "CaloHitCreator/No tail catcher found in GEAR file : " << unknownParameterException.what() << std::endl;
        m_hasTc = false;
    }

    m_mapper = dynamic_cast<const CALICE::AhcMapper*> ( CALICE::MappingProcessor::getMapper(m_settings.m_mappingProcessorName) );
    if (!m_mapper) {
      streamlog_out(ERROR) << "Cannot obtain AhcMapper from MappingProcessor " << m_settings.m_mappingProcessorName
                           <<". Mapper not present or wrong type." << std::endl;
    }
 
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitCreator::~CaloHitCreator()
{
    delete m_pPathLengthCalculator;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode CaloHitCreator::CreateCaloHits(const EVENT::LCEvent *const pLCEvent)
{
    UTIL::CellIDDecoder<CalorimeterHit>::setDefaultEncoding("K-1");

    if (m_hasECal)
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateECalCaloHits(pLCEvent));
    if (m_hasHCal)
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateHCalCaloHits(pLCEvent));
    if (m_hasTc)
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateTcCaloHits(pLCEvent));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode CaloHitCreator::CreateECalCaloHits(const EVENT::LCEvent *const pLCEvent)
{
    for (StringVector::const_iterator iter = m_settings.m_eCalCaloHitCollections.begin(), iterEnd = m_settings.m_eCalCaloHitCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pCaloHitCollection = pLCEvent->getCollection(*iter);
            const int nElements(pCaloHitCollection->getNumberOfElements());

            if (0 == nElements)
                continue;

            static const gear::LayerLayout &layerLayout(marlin::Global::GEAR->getEcalEndcapParameters().getLayerLayout()); 

            UTIL::CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);

//            float xmean = 0.f, ymean =0.f, esum=0.f;

            for (int i = 0; i < nElements; ++i)
            {
                try
                {
                    EVENT::CalorimeterHit *pCaloHit = dynamic_cast<CalorimeterHit*>(pCaloHitCollection->getElementAt(i));

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::ECAL;
                    caloHitParameters.m_isDigital = false;
//                     caloHitParameters.m_layer = cellIdDecoder(pCaloHit)["K-1"]+1;
                    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)["K-1"];
                    caloHitParameters.m_isInOuterSamplingLayer = false;


//                     if (caloHitParameters.m_layer.Get()<=5)
//                     {
//                         const float *pCaloHitPosition(pCaloHit->getPosition());
//                         const pandora::CartesianVector positionVector(pCaloHitPosition[0], pCaloHitPosition[1], pCaloHitPosition[2]);

//                         float hitEnergy = pCaloHit->getEnergy();
//                         esum += hitEnergy;
//                         xmean += positionVector.GetX()*hitEnergy;
//                         ymean += positionVector.GetY()*hitEnergy;
//                     }




                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                    float absorberCorrection(1.);

                    this->GetCaloHitProperties(pCaloHit, layerLayout, caloHitParameters, absorberCorrection, m_settings.m_absorberRadiationLengthECal, m_settings.m_absorberInteractionLengthECal, true);

                    caloHitParameters.m_hadronicEnergy = m_settings.m_eCalToHadGeV * pCaloHit->getEnergy();
                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_eCalToMip * absorberCorrection;

                    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_eCalMipThreshold)
                        continue;

                    caloHitParameters.m_electromagneticEnergy = m_settings.m_eCalToEMGeV * pCaloHit->getEnergy();

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters));
                    m_calorimeterHitVector.push_back(pCaloHit);
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract ecal calo hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (EVENT::Exception &exception)
                {
                    streamlog_out(WARNING) << "Failed to extract ecal calo hit: " << exception.what() << std::endl;
                }

            }
//             if (esum>0)
//             {
//                 xmean /= esum;
//                 ymean /= esum;

//                 std::cout << "ecal x " << xmean << " y " << ymean << std::endl;
//             }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(MESSAGE) << "Failed to extract ecal calo hit collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }
    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode CaloHitCreator::CreateHCalCaloHits(const EVENT::LCEvent *const pLCEvent)
{
    std::cout << "create hcal hits" << std::endl;
    for (StringVector::const_iterator iter = m_settings.m_hCalCaloHitCollections.begin(), iterEnd = m_settings.m_hCalCaloHitCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pCaloHitCollection = pLCEvent->getCollection(*iter);
            const int nElements(pCaloHitCollection->getNumberOfElements());

            const EVENT::LCParameters &collParameters = pCaloHitCollection->getParameters ();
            bool useGearForHCal = collParameters.getNInt("UseGearForHCalCellSize")==1 ? true : false;
            std::cout << "useGearForHCal : " << (useGearForHCal?"yes":"no") << std::endl;
            
//             bool useGearForHCal = false;

            std::cout << "number of hcal hits : " << nElements << std::endl;



            if (0 == nElements)
                continue;

            static const gear::LayerLayout &layerLayout(marlin::Global::GEAR->getHcalEndcapParameters().getLayerLayout());

            UTIL::CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);

            for (int i = 0; i < nElements; ++i)
            {
                try
                {
                    EVENT::CalorimeterHit *pCaloHit = dynamic_cast<CalorimeterHit*>(pCaloHitCollection->getElementAt(i));
                    
                    PandoraApi::CaloHit::Parameters caloHitParameters;

                    caloHitParameters.m_hitType = pandora::HCAL;
                    caloHitParameters.m_isDigital = false;
//                    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)["K-1"]+1;
                    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)["K-1"];
                    caloHitParameters.m_isInOuterSamplingLayer = int(layerLayout.getNLayers()-caloHitParameters.m_layer.Get()) <= m_settings.m_nOuterSamplingLayers;
                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                    float absorberCorrection(1.);

                    this->GetCaloHitProperties(pCaloHit, layerLayout, caloHitParameters, absorberCorrection, m_settings.m_absorberRadiationLengthHCal, m_settings.m_absorberInteractionLengthHCal, useGearForHCal);

                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_hCalToMip * absorberCorrection;

                    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_hCalMipThreshold)
                        continue;

                    caloHitParameters.m_hadronicEnergy = std::min(m_settings.m_hCalToHadGeV * pCaloHit->getEnergy(), m_settings.m_maxHCalHitHadronicEnergy);
                    caloHitParameters.m_electromagneticEnergy = m_settings.m_hCalToEMGeV * pCaloHit->getEnergy();

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters));
                    m_calorimeterHitVector.push_back(pCaloHit);
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract hcal calo hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (EVENT::Exception &exception)
                {
                    streamlog_out(WARNING) << "Failed to extract hcal calo hit: " << exception.what() << std::endl;
                }
            }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(MESSAGE) << "Failed to extract hcal calo hit collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}


//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode CaloHitCreator::CreateTcCaloHits(const EVENT::LCEvent *const pLCEvent)
{
    std::cout << "create tail catcher hits" << std::endl;
    for (StringVector::const_iterator iter = m_settings.m_tcCaloHitCollections.begin(), iterEnd = m_settings.m_tcCaloHitCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pCaloHitCollection = pLCEvent->getCollection(*iter);
            const int nElements(pCaloHitCollection->getNumberOfElements());

//            std::cout << "number of tc hits : " << nElements << std::endl;


            if (0 == nElements)
                continue;

            static const gear::LayerLayout &layerLayout(marlin::Global::GEAR->getYokeEndcapParameters().getLayerLayout());

            UTIL::CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);

            for (int i = 0; i < nElements; ++i)
            {
                try
                {
                    EVENT::CalorimeterHit *pCaloHit = dynamic_cast<CalorimeterHit*>(pCaloHitCollection->getElementAt(i));
                    
                    PandoraApi::CaloHit::Parameters caloHitParameters;


                    if (m_settings.m_makeTcHitsMuonHits)
                        caloHitParameters.m_hitType = pandora::MUON; // you might have to turn of the muon energy corrections, since for the test beam some things are not defined such as the coil
                    else
                        caloHitParameters.m_hitType = pandora::HCAL; 

                    caloHitParameters.m_isDigital = false;
//                    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)["K-1"]+1;
                    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)["K-1"];
                    caloHitParameters.m_isInOuterSamplingLayer = int(layerLayout.getNLayers()-caloHitParameters.m_layer.Get()) <= m_settings.m_nOuterSamplingLayers;

                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                    float absorberCorrection(1.);


                    this->GetCaloHitProperties(pCaloHit, layerLayout, caloHitParameters, absorberCorrection, m_settings.m_absorberRadiationLengthTc, m_settings.m_absorberInteractionLengthTc, true);

                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_tcToMip * absorberCorrection;

                    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_tcMipThreshold)
                        continue;

                    caloHitParameters.m_hadronicEnergy = std::min(m_settings.m_tcToHadGeV * pCaloHit->getEnergy(), m_settings.m_maxTcHitHadronicEnergy);
                    caloHitParameters.m_electromagneticEnergy = m_settings.m_tcToEMGeV * pCaloHit->getEnergy();

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters));
                    m_calorimeterHitVector.push_back(pCaloHit);
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract tail catcher calo hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (EVENT::Exception &exception)
                {
                    streamlog_out(WARNING) << "Failed to extract tail catcher calo hit: " << exception.what() << std::endl;
                }
            }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(MESSAGE) << "Failed to extract tail catcher calo hit collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}



//------------------------------------------------------------------------------------------------------------------------------------------

void CaloHitCreator::GetCommonCaloHitProperties(const EVENT::CalorimeterHit *const pCaloHit, PandoraApi::CaloHit::Parameters &caloHitParameters) const
{
    const float *pCaloHitPosition(pCaloHit->getPosition());
    const pandora::CartesianVector positionVector(pCaloHitPosition[0], pCaloHitPosition[1], pCaloHitPosition[2]);

    caloHitParameters.m_positionVector = positionVector;
    caloHitParameters.m_expectedDirection = positionVector.GetUnitVector();
    caloHitParameters.m_pParentAddress = (void*)pCaloHit;
    caloHitParameters.m_inputEnergy = pCaloHit->getEnergy();
    caloHitParameters.m_time = pCaloHit->getTime();

    // TODO When available, use gear to calculate path length properties.
    // const gear::Vector3D positionIP3D(0.f, 0.f, 0.f);
    // const gear::Vector3D positionVector3D(pCaloHitPosition[0], pCaloHitPosition[1], pCaloHitPosition[2]);
    // radiationLengthsFromIp = marlin::Global::GEAR->getDistanceProperties().getNRadlen(positionIP3D, positionVector3D);
    // interactionLengthsFromIp = marlin::Global::GEAR->getDistanceProperties().getNIntlen(positionIP3D, positionVector3D);

    float radiationLengthsFromIp(0.f), interactionLengthsFromIp(0.f);

    try
    {
        m_pPathLengthCalculator->GetPathLengths(pCaloHit, radiationLengthsFromIp, interactionLengthsFromIp);
    }
    catch (pandora::StatusCodeException &statusCodeException)
    {
        streamlog_out(ERROR) << "Failed to calculate calo hit path lengths from ip: " << statusCodeException.ToString() << std::endl;
        throw statusCodeException;
    }

    try
    {
        caloHitParameters.m_nRadiationLengthsFromIp = radiationLengthsFromIp;
        caloHitParameters.m_nInteractionLengthsFromIp = interactionLengthsFromIp;
    }
    catch (pandora::StatusCodeException &statusCodeException)
    {
        streamlog_out(ERROR) << "Failed to set the radiation or interaction lengths from IP. Values are: radiationLengthsFromIp=" << radiationLengthsFromIp 
                             << "   interactionLengthsFromIp=" << interactionLengthsFromIp << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloHitCreator::GetCaloHitProperties(const EVENT::CalorimeterHit *const pCaloHit, const gear::LayerLayout &layerLayout,
                                          PandoraApi::CaloHit::Parameters &caloHitParameters, float &absorberCorrection, const float absorberRadiationLength, const float absorberInteractionLength, 
                                          bool getCellSizesFromGear) const
{
    caloHitParameters.m_detectorRegion = pandora::ENDCAP;

    int cellID = pCaloHit->getCellID0();
    const int physicalLayer(std::min(static_cast<int>(caloHitParameters.m_layer.Get()), layerLayout.getNLayers() )); //- 1));
    const float layerAbsorberThickness(layerLayout.getAbsorberThickness(physicalLayer));
    caloHitParameters.m_cellThickness = layerLayout.getThickness(physicalLayer); 


    // correction for cell center
    const float zPosCorrection = layerAbsorberThickness/2.f;
    const pandora::CartesianVector positionVector(caloHitParameters.m_positionVector.Get().GetX(),caloHitParameters.m_positionVector.Get().GetY(),caloHitParameters.m_positionVector.Get().GetZ()+zPosCorrection);
    caloHitParameters.m_positionVector = positionVector;


    // get cell dimensions in U and V
    if (getCellSizesFromGear)
    {
        float cellSize0 = layerLayout.getCellSize0(physicalLayer);
        float cellSize1 = layerLayout.getCellSize1(physicalLayer);

//        std::cout << "physical layer : " << physicalLayer << "  thickness : " << caloHitParameters.m_cellThickness.Get() << "  cellSize0 : " << cellSize0 << "   cellSize1 : " << cellSize1 << "  thickness : " << caloHitParameters.m_cellThickness.Get() << "   z pos : " << caloHitParameters.m_positionVector.Get().GetZ() << std::endl;

        caloHitParameters.m_cellSizeU = cellSize1;
        caloHitParameters.m_cellSizeV = cellSize0;
    }
    else
    {
        caloHitParameters.m_cellSizeU = m_mapper->getISizeFromCellID (cellID)*10.0;
        caloHitParameters.m_cellSizeV = m_mapper->getJSizeFromCellID (cellID)*10.0;
    }
    

    caloHitParameters.m_nCellRadiationLengths = absorberRadiationLength * layerAbsorberThickness;
    caloHitParameters.m_nCellInteractionLengths = absorberInteractionLength * layerAbsorberThickness;

    absorberCorrection = 1.;
    for (unsigned int i = 0, iMax = layerLayout.getNLayers(); i < iMax; ++i)
    {
        const float absorberThickness(layerLayout.getAbsorberThickness(i));

        if (absorberThickness <= 0.)
            continue;

        if (layerAbsorberThickness > 0.)
            absorberCorrection = absorberThickness / layerAbsorberThickness;

        break;
    }


    caloHitParameters.m_cellNormalVector = (pCaloHit->getPosition()[2] > 0) ? pandora::CartesianVector(0, 0, 1) :
        pandora::CartesianVector(0, 0, -1);
}


