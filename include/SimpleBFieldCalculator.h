/**
 *  @file   MarlinPandora/include/SimpleBFieldCalculator.h
 * 
 *  @brief  Header file for the simple bfield calculator class.
 * 
 *  $Log: $
 */

#ifndef SIMPLE_BFIELD_CALCULATOR_H
#define SIMPLE_BFIELD_CALCULATOR_H 1

#include "Utilities/BFieldCalculator.h"

/**
 *  @brief  SimpleBFieldCalculator class
 */
class SimpleBFieldCalculator : public pandora::BFieldCalculator
{
public:
    static float        m_bField;                   ///< The bfield in the whole area of the testbeam

private:
    void Initialize(const pandora::GeometryHelper *const pGeometryHelper) {}
    float GetBField(const pandora::CartesianVector &positionVector) const;

};

#endif // #ifndef SIMPLE_BFIELD_CALCULATOR_H
