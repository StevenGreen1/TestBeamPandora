/**
 *  @file   MarlinPandora/src/PathLengthCalculator.cc
 * 
 *  @brief  Implementation of the calo hit creator class.
 * 
 *  $Log: $
 */

#include "marlin/Processor.h"

#include "PathLengthCalculator.h"

#include <limits>

PathLengthCalculator *PathLengthCalculator::GetInstance()
{
    if (!m_instanceFlag)
    {
        m_pPathLengthCalculator = new PathLengthCalculator();
        m_instanceFlag = true;
    }

    return m_pPathLengthCalculator;
}

//------------------------------------------------------------------------------------------------------------------------------------------

PathLengthCalculator::PathLengthCalculator()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

PathLengthCalculator::~PathLengthCalculator()
{
    m_instanceFlag = false;
    m_isGeometryInitialized = false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PathLengthCalculator::GetPathLengths(const EVENT::CalorimeterHit *const pCaloHit, float &nRadiationLengthsFromIp,
    float &nInteractionLengthsFromIp)
{
    if (!m_isGeometryInitialized)
        PathLengthCalculator::InitializeGeometry();

    const pandora::CartesianVector positionVector(pCaloHit->getPosition()[0], pCaloHit->getPosition()[1], pCaloHit->getPosition()[2]);
    const float lineLength(positionVector.GetMagnitude());

    const float eCalPathLength(lineLength * GetFractionInSubDetector(positionVector, m_eCalParameters));
    const float hCalPathLength(lineLength * GetFractionInSubDetector(positionVector, m_hCalParameters));
    const float tcPathLength  (lineLength * GetFractionInSubDetector(positionVector, m_tcParameters  ));

    nRadiationLengthsFromIp = (eCalPathLength * Settings::m_avgRadLengthECal) +
        (hCalPathLength * Settings::m_avgRadLengthHCal) +
        (tcPathLength   * Settings::m_avgRadLengthTc);

    nInteractionLengthsFromIp = (eCalPathLength * Settings::m_avgIntLengthECal) +
        (hCalPathLength * Settings::m_avgIntLengthHCal) +
        (tcPathLength   * Settings::m_avgIntLengthTc);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float PathLengthCalculator::GetFractionInSubDetector(const pandora::CartesianVector &positionVector, const SubDetectorParameters &subDetectorParameters)
{
    const float rCoordinate(std::sqrt(positionVector.GetX() * positionVector.GetX() + positionVector.GetY() * positionVector.GetY()));
    const float zCoordinate(std::fabs(positionVector.GetZ()));

//     if ((rCoordinate < subDetectorParameters.GetInnerRCoordinate()) && (zCoordinate < subDetectorParameters.GetInnerZCoordinate()))
//         return 0.f;

    if (zCoordinate < subDetectorParameters.GetInnerZCoordinate())
        return 0.f;

    // Use r over z ratios to determine where line enters/exits the subdetector
    const float rOverZ((zCoordinate == 0.f) ? std::numeric_limits<float>::max() : rCoordinate / zCoordinate);

//     const float innerROverZ((subDetectorParameters.GetInnerZCoordinate() == 0.f) ? std::numeric_limits<float>::max() :
//         subDetectorParameters.GetInnerRCoordinate() / subDetectorParameters.GetInnerZCoordinate());

    const float outerROverZ((subDetectorParameters.GetOuterZCoordinate() == 0.f) ? std::numeric_limits<float>::max() :
        subDetectorParameters.GetOuterRCoordinate() / subDetectorParameters.GetOuterZCoordinate());

    // Find point at which line enters subdetector
    float innerFraction(0.f);

    innerFraction = GetLengthFraction(positionVector, subDetectorParameters.GetInnerZCoordinate(), subDetectorParameters.GetZNormalVectors());

    if ((innerFraction <= 0.f) || (innerFraction >= 1.f))
        return 0.f;

    // Point at which line exits subdetector
    float outerFraction(0.f);

    if (rOverZ >= outerROverZ)
    {
        outerFraction = GetLengthFraction(positionVector, subDetectorParameters.GetOuterRCoordinate(), subDetectorParameters.GetOuterRNormalVectors());
    }
    else
    {
        outerFraction = GetLengthFraction(positionVector, subDetectorParameters.GetOuterZCoordinate(), subDetectorParameters.GetZNormalVectors());
    }

    if (outerFraction < innerFraction)
        return 0.f;

    return (std::min(1.f, outerFraction) - innerFraction);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float PathLengthCalculator::GetLengthFraction(const pandora::CartesianVector &positionVector, const float closestDistanceToIp,
    const NormalVectorList &normalVectorList)
{
    // Deal with cylindrical case
    if (normalVectorList.empty())
    {
        const float radius(std::sqrt(positionVector.GetX() * positionVector.GetX() + positionVector.GetY() * positionVector.GetY()));
        if (radius<=0)
            return 0.f;
        return closestDistanceToIp / radius;
    }

    // Deal with regular polygon case
    float maxDotProduct(0.f);

    for (NormalVectorList::const_iterator iter = normalVectorList.begin(), iterEnd = normalVectorList.end(); iter != iterEnd; ++iter)
    {
        const float dotProduct(positionVector.GetDotProduct(*iter));

//        std::cout << "pos : " << positionVector.GetX() << ", " << positionVector.GetY() << ", " << positionVector.GetZ() << "  norm " << (*iter).GetX() << ", " << (*iter).GetY() << ", " << (*iter).GetZ() << "   dot prod " << dotProduct << std::endl;

        if (dotProduct > maxDotProduct)
            maxDotProduct = dotProduct;
    }

    if (maxDotProduct <= 0.f)
    {
        streamlog_out(ERROR) << "PathLengthCalculator/GetLengthFraction : maxDotProduct <= 0.f (maxDotProduct=" << maxDotProduct << ") " << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
    }

    return closestDistanceToIp / maxDotProduct;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PathLengthCalculator::InitializeGeometry()
{
    try
    {
        m_eCalParameters.Initialize(pandora::GeometryHelper::GetECalEndCapParameters());
    }
    catch (pandora::StatusCodeException &statusCodeException)
    {
        streamlog_out(MESSAGE) << "PathLengthCalculator: ECal geometry could not be initialized, probably no ECal is present: " << statusCodeException.ToString() << std::endl;
    }

    try
    {
        m_hCalParameters.Initialize(pandora::GeometryHelper::GetHCalEndCapParameters());
        m_isGeometryInitialized = true;
    }
    catch (pandora::StatusCodeException &statusCodeException)
    {
        streamlog_out(ERROR) << "PathLengthCalculator: Failed to initialize HCal: " << statusCodeException.ToString() << std::endl;
        throw statusCodeException;
    }

    try
    {
        m_tcParameters.Initialize(pandora::GeometryHelper::GetMuonEndCapParameters());
        m_isGeometryInitialized = true;
    }
    catch (pandora::StatusCodeException &statusCodeException)
    {
        streamlog_out(ERROR) << "PathLengthCalculator: Failed to initialize tail catcher: " << statusCodeException.ToString() << std::endl;
        throw statusCodeException;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void PathLengthCalculator::SubDetectorParameters::Initialize(const pandora::GeometryHelper::SubDetectorParameters &subDetectorParameters)
{
    m_innerRCoordinate = subDetectorParameters.GetInnerRCoordinate();
    m_outerRCoordinate = subDetectorParameters.GetOuterRCoordinate();
    m_innerZCoordinate = subDetectorParameters.GetInnerZCoordinate();
    m_outerZCoordinate = subDetectorParameters.GetOuterZCoordinate();

    this->GetPolygonNormalVectors(subDetectorParameters.GetInnerSymmetryOrder(), subDetectorParameters.GetInnerPhiCoordinate(), m_innerRNormalVectors);
    this->GetPolygonNormalVectors(subDetectorParameters.GetOuterSymmetryOrder(), subDetectorParameters.GetOuterPhiCoordinate(), m_outerRNormalVectors);
    m_zNormalVectors.push_back(pandora::CartesianVector(0, 0, 1));
//    m_zNormalVectors.push_back(pandora::CartesianVector(0, 0, -1));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PathLengthCalculator::SubDetectorParameters::Initialize_Cylinder(const float innerRCoordinate, const float outerRCoordinate,
    const float innerZCoordinate, const float outerZCoordinate)
{
    m_innerRCoordinate = innerRCoordinate;
    m_outerRCoordinate = outerRCoordinate;
    m_innerZCoordinate = innerZCoordinate;
    m_outerZCoordinate = outerZCoordinate;

    m_zNormalVectors.push_back(pandora::CartesianVector(0, 0, 1));
//    m_zNormalVectors.push_back(pandora::CartesianVector(0, 0, -1));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PathLengthCalculator::SubDetectorParameters::GetPolygonNormalVectors(const unsigned int symmetry, const float phi0, NormalVectorList &normalVectorList) const
{
    for (unsigned int iSymmetry = 0; iSymmetry < symmetry; ++iSymmetry)
    {
        static const float pi(std::acos(-1.));
        const float phi = phi0 + (2. * pi * static_cast<float>(iSymmetry) / static_cast<float>(symmetry));
        normalVectorList.push_back(pandora::CartesianVector(std::sin(phi), std::cos(phi), 0));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PathLengthCalculator::m_instanceFlag = false;
bool PathLengthCalculator::m_isGeometryInitialized = false;
PathLengthCalculator* PathLengthCalculator::m_pPathLengthCalculator = NULL;

PathLengthCalculator::SubDetectorParameters PathLengthCalculator::m_eCalParameters;
PathLengthCalculator::SubDetectorParameters PathLengthCalculator::m_hCalParameters;
PathLengthCalculator::SubDetectorParameters PathLengthCalculator::m_tcParameters;

float PathLengthCalculator::Settings::m_avgRadLengthECal = 0.f;
float PathLengthCalculator::Settings::m_avgRadLengthHCal = 0.f;
float PathLengthCalculator::Settings::m_avgRadLengthTc   = 0.f;

float PathLengthCalculator::Settings::m_avgIntLengthECal = 0.f;
float PathLengthCalculator::Settings::m_avgIntLengthHCal = 0.f;
float PathLengthCalculator::Settings::m_avgIntLengthTc   = 0.f;
