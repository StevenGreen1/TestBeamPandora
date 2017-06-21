/**
 *  @file   MarlinPandora/include/CaloHitCreator.h
 * 
 *  @brief  Header file for the calo hit creator class.
 * 
 *  $Log: $
 */

#ifndef CALO_HIT_CREATOR_H
#define CALO_HIT_CREATOR_H 1

#include "EVENT/CalorimeterHit.h"

#include "gear/LayerLayout.h"

#include "Api/PandoraApi.h"


// CALICE includes
#include "AhcMapper.hh"

#include "MappedContainer.hh"



typedef std::vector<CalorimeterHit *> CalorimeterHitVector;

class PathLengthCalculator;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  CaloHitCreator class
 */
class CaloHitCreator
{
public:

    const CALICE::AhcMapper *m_mapper;

    typedef std::vector<std::string> StringVector;

    /**
     *  @brief  Settings class
     */
    class Settings
    {
    public:
        StringVector    m_eCalCaloHitCollections;               ///< The ecal calorimeter hit collections
        StringVector    m_hCalCaloHitCollections;               ///< The hcal calorimeter hit collections
        StringVector    m_tcCaloHitCollections;                 ///< The tc   calorimeter hit collections

        float           m_absorberRadiationLengthECal;          ///< The absorber radiation length
        float           m_absorberInteractionLengthECal;        ///< The absorber interaction length

        float           m_absorberRadiationLengthHCal;          ///< The absorber radiation length
        float           m_absorberInteractionLengthHCal;        ///< The absorber interaction length

        float           m_absorberRadiationLengthTc;            ///< The absorber radiation length
        float           m_absorberInteractionLengthTc;          ///< The absorber interaction length

        float           m_eCalToMip;                            ///< The calibration from deposited ECal energy to mip
        float           m_hCalToMip;                            ///< The calibration from deposited HCal energy to mip
        float           m_tcToMip;                              ///< The calibration from deposited tail catcher energy to mip
        float           m_eCalMipThreshold;                     ///< Threshold for creating calo hits in the ECal, units mip
        float           m_hCalMipThreshold;                     ///< Threshold for creating calo hits in the HCal, units mip
        float           m_tcMipThreshold;                       ///< Threshold for creating calo hits in the tail catcher, units mip


	bool            m_makeTcHitsMuonHits;                   ///< if true, then mark tail catcher hits as "pandora::MUON", if false, mark them as "pandora::HCAL"

        float           m_eCalToEMGeV;                          ///< The calibration from deposited ECal energy to EM energy
        float           m_eCalToHadGeV;                         ///< The calibration from deposited ECal energy to hadronic energy
        float           m_hCalToEMGeV;                          ///< The calibration from deposited HCal energy to EM energy
        float           m_hCalToHadGeV;                         ///< The calibration from deposited HCal energy to hadronic energy
        float           m_tcToEMGeV;                            ///< The calibration from deposited tail catcher energy to EM energy
        float           m_tcToHadGeV;                           ///< The calibration from deposited tail catcher energy to hadronic energy

        float           m_maxHCalHitHadronicEnergy;             ///< The maximum hadronic energy allowed for a single hcal hit
        float           m_maxTcHitHadronicEnergy;               ///< The maximum hadronic energy allowed for a single tail catcher hit

        int             m_nOuterSamplingLayers;                 ///< Number of layers from edge for hit to be flagged as an outer layer hit
        float           m_layersFromEdgeMaxRearDistance;        ///< Maximum number of layers from candidate outer layer hit to rear of detector
	
	std::string     m_mappingProcessorName;                 ///< Name of the mapping processor
    };

    /**
     *  @brief  Constructor
     * 
     *  @param  settings the creator settings
     */
     CaloHitCreator(const Settings &settings);

    /**
     *  @brief  Destructor
     */
     ~CaloHitCreator();

    /**
     *  @brief  Create calo hits
     * 
     *  @param  pLCEvent the lcio event
     */    
    pandora::StatusCode CreateCaloHits(const LCEvent *const pLCEvent);

    /**
     *  @brief  Get the calorimeter hit vector
     * 
     *  @return The calorimeter hit vector
     */
    static const CalorimeterHitVector &GetCalorimeterHitVector();

    /**
     *  @brief  Reset the calo hit creator
     */
    void Reset();

private:
    /**
     *  @brief  Create ecal calo hits
     * 
     *  @param  pLCEvent the lcio event
     */
    pandora::StatusCode CreateECalCaloHits(const EVENT::LCEvent *const pLCEvent);

    /**
     *  @brief  Create hcal calo hits
     * 
     *  @param  pLCEvent the lcio event
     */
    pandora::StatusCode CreateHCalCaloHits(const EVENT::LCEvent *const pLCEvent);

    /**
     *  @brief  Create tail catcher calo hits
     * 
     *  @param  pLCEvent the lcio event
     */
    pandora::StatusCode CreateTcCaloHits(const EVENT::LCEvent *const pLCEvent);

    /**
     *  @brief  Get common calo hit properties: position, parent address, input energy and time
     * 
     *  @param  pCaloHit the lcio calorimeter hit
     *  @param  caloHitParameters the calo hit parameters to populate
     */
    void GetCommonCaloHitProperties(const EVENT::CalorimeterHit *const pCaloHit, PandoraApi::CaloHit::Parameters &caloHitParameters) const;

    /**
     *  @brief  Get calo hit properties: cell size, absorber radiation and interaction lengths, normal vector
     * 
     *  @param  pCaloHit the lcio calorimeter hit
     *  @param  layerLayout the gear end cap layer layout
     *  @param  caloHitParameters the calo hit parameters to populate
     *  @param  absorberCorrection to receive the absorber thickness correction for the mip equivalent energy
     *  @param  getCellSizesFromGear take cell sizes from gear file
     */
    void GetCaloHitProperties(const EVENT::CalorimeterHit *const pCaloHit, const gear::LayerLayout &layerLayout,
			      PandoraApi::CaloHit::Parameters &caloHitParameters, float &absorberCorrection, const float absorberRadiationLength, 
			      const float absorberInteractionLength, bool getCellSizesFromGear=false) const;

    const Settings                      m_settings;                         ///< The calo hit creator settings

    const pandora::Pandora             *m_pPandora;                         ///< Address of the pandora object to create calo hits
    const PathLengthCalculator         *m_pPathLengthCalculator;            ///< Address of the path length calculator


    float                         m_eCalZ;                            ///< ECal z coordinate
    float                         m_eCalR;                            ///< ECal r coordinate
    float                         m_eCalPhi0;                         ///< ECal phi0 coordinate
    unsigned int                  m_eCalSymmetry;                     ///< ECal symmetry order

    float                         m_hCalZ;                            ///< HCal z coordinate
    float                         m_hCalR;                            ///< HCal r coordinate
    float                         m_hCalPhi0;                         ///< HCal phi0 coordinate
    unsigned int                  m_hCalSymmetry;                     ///< HCal symmetry order

    float                         m_tcZ;                              ///< HCal z coordinate
    float                         m_tcR;                              ///< HCal r coordinate
    float                         m_tcPhi0;                           ///< HCal phi0 coordinate
    unsigned int                  m_tcSymmetry  ;                     ///< HCal symmetry order

    bool                          m_hasECal;                          ///< There is an ECal in the detector set up
    bool                          m_hasHCal;                          ///< There is an HCal in the detector set up
    bool                          m_hasTc;                            ///< There is an tail catcher in the detector set up

    float                               m_eCalLayerThickness;               ///< ECal layer thickness
    float                               m_hCalLayerThickness;               ///< HCal layer thickness
    float                               m_tcLayerThickness;                 ///< tail catcher layer thickness

    static CalorimeterHitVector         m_calorimeterHitVector;             ///< The calorimeter hit vector
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline const CalorimeterHitVector &CaloHitCreator::GetCalorimeterHitVector()
{
    return m_calorimeterHitVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void CaloHitCreator::Reset()
{
    m_calorimeterHitVector.clear();
}

#endif // #ifndef CALO_HIT_CREATOR_H
