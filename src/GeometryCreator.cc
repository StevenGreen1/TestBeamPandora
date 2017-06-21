/**
 *  @file   MarlinPandora/src/GeometryCreator.cc
 * 
 *  @brief  Implementation of the geometry creator class.
 * 
 *  $Log: $
 */

#include "marlin/Global.h"
#include "marlin/Processor.h"

#include "gear/GEAR.h"
#include "gear/GearParameters.h"
#include "gear/CalorimeterParameters.h"
#include "gear/PadRowLayout2D.h"
#include "gear/LayerLayout.h"

#include "GeometryCreator.h"
#include "CalicePandoraProcessor.h"

GeometryCreator::GeometryCreator(const Settings &settings) :
    m_settings(settings),
    m_pPandora(CalicePandoraProcessor::GetPandora())
{
    streamlog_out(MESSAGE) << "GeometryCreator - c'tor" << std::endl; // BBB0
}

//------------------------------------------------------------------------------------------------------------------------------------------

GeometryCreator::~GeometryCreator()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode GeometryCreator::CreateGeometry() const
{
    streamlog_out(MESSAGE) << "GeometryCreator - Create Geometry" << std::endl; // BBB0

    PandoraApi::Geometry::Parameters geometryParameters;

    // Initialize settings to gear defaults
    try
    {
        const gear::CalorimeterParameters &eCalEndcapParameters = marlin::Global::GEAR->getEcalEndcapParameters();
        SetDefaultSubDetectorParameters(eCalEndcapParameters, geometryParameters.m_eCalEndCapParameters, m_settings.m_absorberRadiationLengthECal, m_settings.m_absorberInteractionLengthECal );
    }
    catch (gear::Exception &exception)
    {
        streamlog_out(MESSAGE) << "CreateGeometry/No ECal Endcap parameters found in GEAR file: " << exception.what() << std::endl;
    }

    try
    {
        const gear::CalorimeterParameters &hCalEndcapParameters = marlin::Global::GEAR->getHcalEndcapParameters();
        SetDefaultSubDetectorParameters(hCalEndcapParameters, geometryParameters.m_hCalEndCapParameters, m_settings.m_absorberRadiationLengthHCal, m_settings.m_absorberInteractionLengthHCal);
    }
    catch (gear::Exception &exception)
    {
        streamlog_out(MESSAGE) << "CreateGeometry/No HCal Endcap parameters found in GEAR file: " << exception.what() << std::endl;
    }

    try
    {
        const gear::CalorimeterParameters &YokeEndcapParameters = marlin::Global::GEAR->getYokeEndcapParameters();
        SetDefaultSubDetectorParameters(YokeEndcapParameters, geometryParameters.m_muonEndCapParameters, m_settings.m_absorberRadiationLengthTc, m_settings.m_absorberInteractionLengthTc);
    }
    catch (gear::Exception &exception)
    {
        streamlog_out(MESSAGE) << "CreateGeometry/No tail catcher parameters found in GEAR file: " << exception.what() << std::endl;
    }


    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::Create(*m_pPandora, geometryParameters));

   return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GeometryCreator::SetDefaultSubDetectorParameters(const gear::CalorimeterParameters &inputParameters,
    PandoraApi::GeometryParameters::SubDetectorParameters &subDetectorParameters, const float absorberRadiationLength, const float absorberInteractionLength) const
{
    const gear::LayerLayout &layerLayout = inputParameters.getLayerLayout();

    subDetectorParameters.m_innerRCoordinate    = inputParameters.getExtent()[0];
    subDetectorParameters.m_innerZCoordinate    = inputParameters.getExtent()[2];
    subDetectorParameters.m_innerPhiCoordinate  = inputParameters.getPhi0();
    subDetectorParameters.m_innerSymmetryOrder  = inputParameters.getSymmetryOrder();
    subDetectorParameters.m_outerRCoordinate    = inputParameters.getExtent()[1];
    subDetectorParameters.m_outerZCoordinate    = inputParameters.getExtent()[3];
    subDetectorParameters.m_outerPhiCoordinate  = inputParameters.getPhi0();
    subDetectorParameters.m_outerSymmetryOrder  = inputParameters.getSymmetryOrder();
    subDetectorParameters.m_isMirroredInZ       = false;
    subDetectorParameters.m_nLayers             = layerLayout.getNLayers();

    for (int i = 0; i < layerLayout.getNLayers(); ++i)
    {
        PandoraApi::Geometry::Parameters::LayerParameters layerParameters;
        layerParameters.m_closestDistanceToIp   = layerLayout.getDistance(i) + (0.5 * (layerLayout.getThickness(i) + layerLayout.getAbsorberThickness(i)));
        layerParameters.m_nRadiationLengths     = absorberRadiationLength * layerLayout.getAbsorberThickness(i);
        layerParameters.m_nInteractionLengths   = absorberInteractionLength * layerLayout.getAbsorberThickness(i);
        subDetectorParameters.m_layerParametersList.push_back(layerParameters);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GeometryCreator::SetAdditionalSubDetectorParameters(PandoraApi::GeometryParameters &geometryParameters) const
{
// add for each additional subdetector a section like the following example (which is commented out)

//     try
//     {
//         PandoraApi::Geometry::Parameters::SubDetectorParameters eCalPlugParameters;
//         const gear::CalorimeterParameters &eCalPlugInputParameters = marlin::Global::GEAR->getEcalPlugParameters();
//         SetDefaultSubDetectorParameters(eCalPlugInputParameters, eCalPlugParameters);
//         geometryParameters.m_additionalSubDetectors["ECalPlug"] = eCalPlugParameters;
//     }
//     catch (gear::Exception &exception)
//     {
//         streamlog_out(WARNING) << "Marlin pandora geometry creator: " << exception.what() << std::endl;
//     }

}


