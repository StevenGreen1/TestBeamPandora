/**
 *  @file   PandoraPFANew/src/Utilities/CalicePseudoLayerCalculator.cc
 * 
 *  @brief  Implementation of the fine granularity pseudo layer calculator class.
 * 
 *  $Log: $
 */

#include "CalicePseudoLayerCalculator.h"


#include "marlin/Global.h"
#include "marlin/Processor.h"

#include "gear/GearParameters.h"
#include "gear/CalorimeterParameters.h"
#include "gear/GEAR.h"
#include "gear/LayerLayout.h"


#include <algorithm>

using namespace pandora;

CalicePseudoLayerCalculator::CalicePseudoLayerCalculator() :
    PseudoLayerCalculator()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalicePseudoLayerCalculator::InitializeGeometry()
{
    float layerStart = 0.0;
    int pseudoLayerNumber = 0;
    try
    {
        const float eCalZ = marlin::Global::GEAR->getEcalEndcapParameters().getExtent()[2];
        const gear::LayerLayout &eCalLayerLayout(marlin::Global::GEAR->getEcalEndcapParameters().getLayerLayout()); 

        layerStart = eCalZ;
        m_layerPositions.insert (std::make_pair(layerStart,pseudoLayerNumber));
        ++pseudoLayerNumber;
        for (int layer = 0, layerEnd = eCalLayerLayout.getNLayers(); layer < layerEnd; ++layer)
        {
            const float eCalLayerThickness = eCalLayerLayout.getThickness (layer);
            layerStart += eCalLayerThickness;
            m_layerPositions.insert (std::make_pair(layerStart,pseudoLayerNumber));
            ++pseudoLayerNumber;
        }
    }
    catch (gear::UnknownParameterException &unknownParameterException)
    {
        streamlog_out(MESSAGE) << "No ECal found in GEAR file : " << unknownParameterException.what() << std::endl;
    }
    

    try
    {
        const float hCalZ = marlin::Global::GEAR->getHcalEndcapParameters().getExtent()[2];
        const gear::LayerLayout &hCalLayerLayout(marlin::Global::GEAR->getHcalEndcapParameters().getLayerLayout()); 

        layerStart = hCalZ;
        m_layerPositions.insert (std::make_pair(layerStart,pseudoLayerNumber));
        ++pseudoLayerNumber;
        for (int layer = 0, layerEnd = hCalLayerLayout.getNLayers(); layer < layerEnd; ++layer)
        {
            const float hCalLayerThickness = hCalLayerLayout.getThickness (layer);
            layerStart += hCalLayerThickness;
            m_layerPositions.insert (std::make_pair(layerStart,pseudoLayerNumber));
            ++pseudoLayerNumber;
        }
    }
    catch (gear::UnknownParameterException &unknownParameterException)
    {
        streamlog_out(MESSAGE) << "No HCal found in GEAR file : " << unknownParameterException.what() << std::endl;
    }

    try
    {
        const float tcZ = marlin::Global::GEAR->getYokeEndcapParameters().getExtent()[2];
        const gear::LayerLayout &tcLayerLayout(marlin::Global::GEAR->getYokeEndcapParameters().getLayerLayout()); 

        layerStart = tcZ;
        m_layerPositions.insert (std::make_pair(layerStart,pseudoLayerNumber));
        ++pseudoLayerNumber;
        for (int layer = 0, layerEnd = tcLayerLayout.getNLayers(); layer < layerEnd; ++layer)
        {
            const float tcLayerThickness = tcLayerLayout.getThickness (layer);
            layerStart += tcLayerThickness;
            m_layerPositions.insert (std::make_pair(layerStart,pseudoLayerNumber));
            ++pseudoLayerNumber;
        }
    }
    catch (gear::UnknownParameterException &unknownParameterException)
    {
        streamlog_out(MESSAGE) << "No tail catcher found in GEAR file : " << unknownParameterException.what() << std::endl;
    }


    for (LayerPositionMap::iterator itPos = m_layerPositions.begin(), itPosEnd = m_layerPositions.end(); itPos != itPosEnd; ++itPos)
    {
        float start = itPos->first;
        int   number = itPos->second;
        std::cout << "pseudo layer " << number << " starts at " << start << std::endl;
    }

}

//------------------------------------------------------------------------------------------------------------------------------------------

PseudoLayer CalicePseudoLayerCalculator::GetPseudoLayer(const CartesianVector &positionVector) const
{
    const float zCoordinate(std::fabs(positionVector.GetZ()));

    LayerPositionMap::const_iterator itLayerPosition = --(m_layerPositions.lower_bound (zCoordinate));
    const PseudoLayer pseudoLayer = (*itLayerPosition).second;
//    std::cout << "get pseudo layer for " << zCoordinate << " : layer " << pseudoLayer << std::endl;

    // Reserve a pseudo layer for track projections, etc.
    return (1 + pseudoLayer);
}

