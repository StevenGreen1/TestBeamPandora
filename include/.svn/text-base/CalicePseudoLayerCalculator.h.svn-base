/**
 *  @file   PandoraPFANew/include/Utilities/CalicePseudoLayerCalculator.h
 * 
 *  @brief  Header file for the fine granularity pseudo layer calculator class.
 * 
 *  $Log: $
 */
#ifndef CALICE_PSEUDO_LAYER_CALCULATOR_H
#define CALICE_PSEUDO_LAYER_CALCULATOR_H 1

#include "Helpers/GeometryHelper.h"

#include "Utilities/PseudoLayerCalculator.h"

/**
 *  @brief  CalicePseudoLayerCalculator class
 */
class CalicePseudoLayerCalculator : public pandora::PseudoLayerCalculator
{
public:
    /**
     *  @brief  Default constructor
     */
    CalicePseudoLayerCalculator();

private:
    void InitializeGeometry();
    pandora::PseudoLayer GetPseudoLayer(const pandora::CartesianVector &positionVector) const;
    pandora::PseudoLayer GetPseudoLayerAtIp() const;

    typedef std::map<float,int> LayerPositionMap;

    LayerPositionMap m_layerPositions; ///< list with the positions of the layers
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::PseudoLayer CalicePseudoLayerCalculator::GetPseudoLayerAtIp() const
{
    static const pandora::PseudoLayer pseudoLayerAtIp(this->GetPseudoLayer(pandora::CartesianVector(0.f, 0.f, 0.f)));
    return pseudoLayerAtIp;
}

#endif // #ifndef CALICE_PSEUDO_LAYER_CALCULATOR_H
