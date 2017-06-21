/**
 *  @file   MarlinPandora/include/TrackCreator.h
 * 
 *  @brief  Header file for the track creator class.
 * 
 *  $Log: $
 */

#ifndef TRACK_CREATOR_H
#define TRACK_CREATOR_H 1

#include "EVENT/LCEvent.h"
#include "EVENT/Track.h"

#include "TBTrackUtil/TrackProjection.hh"

#include "Api/PandoraApi.h"
#include "Objects/Helix.h"

typedef std::vector<Track *> TrackVector;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TrackCreator class
 */
class TrackCreator
{
public:
    typedef std::vector<double> DoubleVector;
    typedef std::vector<std::string> StringVector;

    /**
     *  @brief  Settings class
     */
    class Settings
    {
    public:
	StringVector    m_trackCollectionsX;                     ///< The reconstructed track collections for X
	StringVector    m_trackCollectionsY;                     ///< The reconstructed track collections for Y

/* 	std::string    m_trackCollectionX;                      ///< The reconstructed track collection for X */
/* 	std::string    m_trackCollectionY;                      ///< The reconstructed track collection for Y */

	float           m_trackStartInZ;                        ///< z-coordinate starting point of the track 
	float           m_trackEndInZ;                          ///< z-coordinate end point of the track 

	float           m_momentumMagnitude;                    ///< momentum magnitude of the track

//	StringVector    m_trackCollections;                     ///< track collection to read in
    };


    /**
     *  @brief  CaliceTrackInfo class
     */
    class CaliceTrackInfo
    {
    public:
	TBTrack::TrackProjection*     m_trackProjectionX;                      ///< The reconstructed track projection for X
	TBTrack::TrackProjection*     m_trackProjectionY;                      ///< The reconstructed track projection for Y
    };

    /**
     *  @brief  Constructor
     * 
     *  @param  settings the creator settings
     */
     TrackCreator(const Settings &settings);

    /**
     *  @brief  Destructor
     */
     ~TrackCreator();

    /**
     *  @brief  Create associations between tracks, V0s, kinks, etc
     * 
     *  @param  pLCEvent the lcio event
     */
    pandora::StatusCode CreateTrackAssociations(const EVENT::LCEvent *const pLCEvent);

    /**
     *  @brief  Create tracks, insert user code here
     * 
     *  @param  pLCEvent the lcio event
     */
    pandora::StatusCode CreateTracks(const EVENT::LCEvent *const pLCEvent);
//    pandora::StatusCode CreateTracks(const EVENT::LCEvent *const pLCEvent) const;

    /**
     *  @brief  Reset the track creator
     */
    void Reset();

private:
    /**
     *  @brief  Project helix to a given plane in z
     *
     *  @param  pHelix helix fit to be projected to the calorimeter surface
     *  @param  signPz sign w.r.t. increasing z direction
     *  @param  zPlane z-position of the plane
     *  @param  radius outer radius of the polygon of the calorimeter
     *  @param  phi0 starting angle of the polygon
     *  @param  symmetry of the polygon
     *  @param  trackParameters the track parameters
     */
    void GetZPlaneProjection(const pandora::Helix *const pHelix, const int signPz, const float zPlane, const float radius,
			     const float phi0, const float symmetry, PandoraApi::Track::Parameters &trackParameters) const;

    /**
     *  @brief  Get a track projection from a collection given by the collectionName
     *
     *  @param  pLCEvent the LCEvent
     *  @param  pTrackProjection a Calice track projection is returned
     *  @param  collectionName name of the collection to be taken to compute the track projection (one collection for X and one for Y)
     *  @param  collection returns the LCCollection which corresponds to the collectionName
     */
    void GetTBTrackProjection(const EVENT::LCEvent *const pLCEvent, TBTrack::TrackProjection*& pTrackProjection,
			      const std::string collectionName, LCCollection*& collection) const;


/*     /\** */
/*      *  @brief  Perform helix fits to calculate track parameters: momentum at dca, start and end track states */
/*      *  */
/*      *  @param  pTrack the lcio track */
/*      *  @param  trackParameters the track parameters */
/*      *\/ */
/*     void FitTrackHelices(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const; */


    const Settings          m_settings;                     ///< The track creator settings
    const pandora::Pandora *m_pPandora;                     ///< Address of the pandora object to create tracks and track relationships

    const float             m_bField;                       ///< The bfield

    typedef std::set<CaliceTrackInfo*> CaliceTrackList;
    CaliceTrackList         m_caliceTrackList;              ///< calice track list which stores the track projections

    static TrackVector      m_trackVector;                  ///< The track vector
    
    float                   m_calorimeterFaceZ;             ///< first calorimeter z position
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackCreator::Reset()
{
    m_trackVector.clear();
}



#endif // #ifndef TRACK_CREATOR_H
