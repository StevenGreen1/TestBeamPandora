/**
 *  @file   MarlinPandora/include/GeometryCreator.h
 * 
 *  @brief  Header file for the geometry creator class.
 * 
 *  $Log: $
 */

#ifndef GEOMETRY_CREATOR_H
#define GEOMETRY_CREATOR_H 1

#include "gear/CalorimeterParameters.h"

#include "Api/PandoraApi.h"

/**
 *  @brief  GeometryCreator class
 */
class GeometryCreator
{
public:
    /**
     *  @brief  Settings class
     */
    class Settings
    {
    public:
        float           m_absorberRadiationLengthECal;              ///< The absorber radiation length
        float           m_absorberInteractionLengthECal;            ///< The absorber interaction length

        float           m_absorberRadiationLengthHCal;              ///< The absorber radiation length
        float           m_absorberInteractionLengthHCal;            ///< The absorber interaction length

        float           m_absorberRadiationLengthTc;                ///< The absorber radiation length
        float           m_absorberInteractionLengthTc;              ///< The absorber interaction length

        int             m_eCalSymmetryOrder;         ///< ECal symmetry order (missing from ILD gear files)
        float           m_eCalPhiCoordinate;         ///< ECal phi coordinate (missing from ILD gear files)

        int             m_hCalSymmetryOrder;         ///< HCal symmetry order (missing from ILD gear files)
        float           m_hCalPhiCoordinate;         ///< HCal phi coordinate (missing from ILD gear files)

        int             m_tcSymmetryOrder;           ///< tail catcher symmetry order (missing from ILD gear files)
        float           m_tcPhiCoordinate;           ///< tail catcher phi coordinate (missing from ILD gear files)
    };

    /**
     *  @brief  Constructor
     * 
     *  @param  settings the creator settings
     */
     GeometryCreator(const Settings &settings);

    /**
     *  @brief  Destructor
     */
     ~GeometryCreator();

    /**
     *  @brief  Create geometry
     */
    pandora::StatusCode CreateGeometry() const;

private:
    /**
     *  @brief  Set sub detector parameters to their gear default values
     * 
     *  @param  inputParameters input parameters, from gear
     *  @param  subDetectorParameters the sub detector parameters
     */
    void SetDefaultSubDetectorParameters(const gear::CalorimeterParameters &inputParameters,
        PandoraApi::GeometryParameters::SubDetectorParameters &subDetectorParameters, const float absorberRadiationLength, const float absorberInteractionLength) const;

    /**
     *  @brief  Set additional sub detector parameters
     * 
     *  @param  geometryParameters the pandora geometry parameters
     */
    void SetAdditionalSubDetectorParameters(PandoraApi::GeometryParameters &geometryParameters) const;

    const Settings          m_settings;                     ///< The geometry creator settings
    const pandora::Pandora *m_pPandora;                     ///< Address of the pandora object to create the geometry
};

#endif // #ifndef GEOMETRY_CREATOR_H
