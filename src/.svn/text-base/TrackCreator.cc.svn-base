/**
 *  @file   MarlinPandora/src/TrackCreator.cc
 * 
 *  @brief  Implementation of the track creator class.
 * 
 *  $Log: $
 */

#include "marlin/Global.h"
#include "marlin/Processor.h"

#include "EVENT/LCCollection.h"
#include "EVENT/LCGenericObject.h"

#include "gear/CalorimeterParameters.h"

#include "CalicePandoraProcessor.h"
#include "TrackCreator.h"
#include "Pandora/PdgTable.h"

#include <algorithm>
#include <cmath>
#include <limits>



TrackVector TrackCreator::m_trackVector;


//------------------------------------------------------------------------------------------------------------------------------------------

TrackCreator::TrackCreator(const Settings &settings) :
    m_settings(settings),
    m_pPandora(CalicePandoraProcessor::GetPandora()),
    m_bField(0.0001)
{
    bool hasFoundAnything = false;
    m_calorimeterFaceZ = std::numeric_limits<float>::max();

    try
    {
        m_calorimeterFaceZ = std::min(m_calorimeterFaceZ,float(marlin::Global::GEAR->getEcalEndcapParameters().getExtent()[2]));
        hasFoundAnything = true;
    }
    catch (gear::UnknownParameterException &unknownParameterException)
    {
        streamlog_out(MESSAGE) << "No ECal found in GEAR file : " << unknownParameterException.what() << std::endl;
    }
    

    try
    {
        m_calorimeterFaceZ = std::min(m_calorimeterFaceZ,float(marlin::Global::GEAR->getHcalEndcapParameters().getExtent()[2]));
        hasFoundAnything = true;
    }
    catch (gear::UnknownParameterException &unknownParameterException)
    {
        streamlog_out(MESSAGE) << "No HCal found in GEAR file : " << unknownParameterException.what() << std::endl;
    }

    try
    {
        m_calorimeterFaceZ = std::min(m_calorimeterFaceZ,float(marlin::Global::GEAR->getYokeEndcapParameters().getExtent()[2]));
        hasFoundAnything = true;
    }
    catch (gear::UnknownParameterException &unknownParameterException)
    {
        streamlog_out(MESSAGE) << "No tail catcher found in GEAR file : " << unknownParameterException.what() << std::endl;
    }

    if (!hasFoundAnything)
        streamlog_out(ERROR) << "TrackCreator/no calorimeter has been found, calorimeter front face could not be determined" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackCreator::~TrackCreator()
{
}



//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode TrackCreator::CreateTracks(const EVENT::LCEvent *const pLCEvent) 
{
    int sizeCollectionX=m_settings.m_trackCollectionsX.size();
    int sizeCollectionY=m_settings.m_trackCollectionsY.size();

    if (sizeCollectionX!=sizeCollectionY)
        return pandora::STATUS_CODE_FAILURE;

    
    for (StringVector::const_iterator iterX = m_settings.m_trackCollectionsX.begin(), iterEndX = m_settings.m_trackCollectionsX.end(), 
                                      iterY = m_settings.m_trackCollectionsY.begin(), iterEndY = m_settings.m_trackCollectionsY.end();
                                      (iterX != iterEndX) && (iterY != iterEndY); )
    {
        std::string trackCollectionNameX=(*iterX);
        std::string trackCollectionNameY=(*iterY);


        TBTrack::TrackProjection* pTrackProjectionX = NULL;
        TBTrack::TrackProjection* pTrackProjectionY = NULL;

        LCCollection* collectionX = NULL;
        LCCollection* collectionY = NULL;

        this->GetTBTrackProjection(pLCEvent, pTrackProjectionX, trackCollectionNameX, collectionX);
        this->GetTBTrackProjection(pLCEvent, pTrackProjectionY, trackCollectionNameY, collectionY);

        if (!pTrackProjectionX || !pTrackProjectionY)
            return pandora::STATUS_CODE_SUCCESS;

        PandoraApi::Track::Parameters trackParameters;
    
        trackParameters.m_d0 = 0.f; // TODO: think if 0.f is a sensible number, or if there is a better one
        trackParameters.m_z0 = 0.f;

        CaliceTrackInfo* caliceTrackInfo = new CaliceTrackInfo();
        caliceTrackInfo->m_trackProjectionX = pTrackProjectionX;
        caliceTrackInfo->m_trackProjectionY = pTrackProjectionY;

        if (!m_caliceTrackList.insert(caliceTrackInfo).second)
            return pandora::STATUS_CODE_FAILURE;

        trackParameters.m_pParentAddress = caliceTrackInfo;



        // By default, assume tracks are charged pions
        trackParameters.m_particleId = pandora::PI_PLUS;
        trackParameters.m_mass = pandora::PdgTable::GetParticleMass(pandora::PI_PLUS);

        trackParameters.m_charge = 1;

        pandora::CartesianVector startPosition, startMomentum, calorimeterPosition;
        pandora::CartesianVector endPosition, endMomentum;

        startPosition.SetValues(pTrackProjectionX->intercept(m_settings.m_trackStartInZ,0),pTrackProjectionY->intercept(m_settings.m_trackStartInZ,0),m_settings.m_trackStartInZ);
        endPosition  .SetValues(pTrackProjectionX->intercept(m_settings.m_trackEndInZ,0),  pTrackProjectionY->intercept(m_settings.m_trackEndInZ,0),  m_settings.m_trackEndInZ);
        calorimeterPosition.SetValues(pTrackProjectionX->intercept(m_calorimeterFaceZ,0),  pTrackProjectionY->intercept(m_calorimeterFaceZ,0),        m_calorimeterFaceZ);


        pandora::CartesianVector momentum = endPosition-startPosition;
        momentum = momentum.GetUnitVector();
        momentum *= m_settings.m_momentumMagnitude;

        trackParameters.m_momentumAtDca = momentum;

        trackParameters.m_trackStateAtStart       = pandora::TrackState(startPosition, momentum);
        trackParameters.m_trackStateAtEnd         = pandora::TrackState(endPosition, momentum);
        trackParameters.m_trackStateAtCalorimeter = pandora::TrackState(calorimeterPosition, momentum);

        //     std::cout << "track at calo: " << calorimeterPosition.GetX() << "  " << calorimeterPosition.GetY() << "  " << calorimeterPosition.GetZ() << std::endl;

        trackParameters.m_timeAtCalorimeter     = 0.f;
        trackParameters.m_isProjectedToEndCap   = true;

        trackParameters.m_reachesCalorimeter    = true;
        trackParameters.m_canFormPfo            = true;
        trackParameters.m_canFormClusterlessPfo = false;


        try
        {
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Track::Create(*m_pPandora, trackParameters));
        }
        catch (pandora::StatusCodeException &statusCodeException)
        {
            streamlog_out(ERROR) << "Failed to extract a track: " << statusCodeException.ToString() << std::endl;
        }

        ++iterX;
        ++iterY;
    }

    return pandora::STATUS_CODE_SUCCESS;
}


//------------------------------------------------------------------------------------------------------------------------------------------

void TrackCreator::GetTBTrackProjection(const EVENT::LCEvent *const pLCEvent, TBTrack::TrackProjection*& pTrackProjection, const std::string collectionName, LCCollection*& collection) const
{
    collection = NULL;
    /* get best track projection */
    try
    {
        collection = pLCEvent->getCollection( collectionName );
    }
    catch ( lcio::DataNotAvailableException &e ) 
    {
        streamlog_out( DEBUG ) << "TrackCreator/GetTBTrackProjection/Missing collection " << collectionName << std::endl;
    }   

    int numberElements = collection->getNumberOfElements();
    for (int iElement = 0; iElement < numberElements; ++iElement)
    {
        TBTrack::TrackProjection *pCurrentTrackProjection = 
            new TBTrack::TrackProjection(dynamic_cast< EVENT::LCGenericObject* >(collection->getElementAt(iElement)) );
	  
        if (pTrackProjection == NULL)
            pTrackProjection = pCurrentTrackProjection;
        else
        {
            if (pCurrentTrackProjection->probability() > pTrackProjection->probability())
            {
                delete pTrackProjection;
                pTrackProjection = pCurrentTrackProjection;
            }
            else
                delete pCurrentTrackProjection;
        }
    }
}







// //------------------------------------------------------------------------------------------------------------------------------------------

// pandora::StatusCode TrackCreator::CreateTracks(const EVENT::LCEvent *const pLCEvent) const
// {
//     for (StringVector::const_iterator iter = m_settings.m_trackCollections.begin(), iterEnd = m_settings.m_trackCollections.end();
//         iter != iterEnd; ++iter)
//     {
//         try
//         {
//             const EVENT::LCCollection *pTrackCollection = pLCEvent->getCollection(*iter);

//             for (int i = 0, iMax = pTrackCollection->getNumberOfElements(); i < iMax; ++i)
//             {
//                 try
//                 {
//                     EVENT::Track *pTrack = dynamic_cast<Track*>(pTrackCollection->getElementAt(i));

//                     // Proceed to create the pandora track
//                     PandoraApi::Track::Parameters trackParameters;
//                     trackParameters.m_d0 = pTrack->getD0();
//                     trackParameters.m_z0 = pTrack->getZ0();
//                     trackParameters.m_pParentAddress = pTrack;

//                     // By default, assume tracks are charged pions
//                     const float signedCurvature(pTrack->getOmega());
//                     trackParameters.m_particleId = (signedCurvature > 0) ? pandora::PI_PLUS : pandora::PI_MINUS;
//                     trackParameters.m_mass = pandora::PdgTable::GetParticleMass(pandora::PI_PLUS);

//                     if (0.f != signedCurvature)
//                         trackParameters.m_charge = static_cast<int>(signedCurvature / std::fabs(signedCurvature));

//                     this->FitTrackHelices(pTrack, trackParameters);

//                     PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Track::Create(*m_pPandora, trackParameters));
//                     m_trackVector.push_back(pTrack);
//                 }
//                 catch (pandora::StatusCodeException &statusCodeException)
//                 {
//                     streamlog_out(ERROR) << "Failed to extract a track: " << statusCodeException.ToString() << std::endl;
//                 }
//             }
//         }
//         catch (EVENT::Exception &exception)
//         {
//             streamlog_out(WARNING) << "Failed to extract track collection: " << *iter << ", " << exception.what() << std::endl;
//         }
//     }

//     return pandora::STATUS_CODE_SUCCESS;
// }

// //------------------------------------------------------------------------------------------------------------------------------------------

// void TrackCreator::FitTrackHelices(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
// {
//     pandora::Helix *pHelixFit = new pandora::Helix(pTrack->getPhi(), pTrack->getD0(), pTrack->getZ0(), pTrack->getOmega()+0.0001, pTrack->getTanLambda(), m_bField);
//     trackParameters.m_momentumAtDca = pHelixFit->GetMomentum();

//     pandora::CartesianVector startPosition, startMomentum;
//     PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, pHelixFit->GetPointInZ(m_settings.m_trackStartInZ, pHelixFit->GetReferencePoint(), startPosition));
//     startMomentum = pHelixFit->GetExtrapolatedMomentum(startPosition);
//     trackParameters.m_trackStateAtStart = pandora::TrackState(startPosition, startMomentum);

//     pandora::CartesianVector endPosition, endMomentum;
//     PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, pHelixFit->GetPointInZ(m_settings.m_trackEndInZ, pHelixFit->GetReferencePoint(), endPosition));
//     endMomentum = pHelixFit->GetExtrapolatedMomentum(endPosition);
//     trackParameters.m_trackStateAtEnd = pandora::TrackState(endPosition, endMomentum);

//     pandora::CartesianVector calorimeterPosition, calorimeterMomentum;
//     PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, pHelixFit->GetPointInZ(m_calorimeterFaceZ, pHelixFit->GetReferencePoint(), calorimeterPosition));
//     calorimeterMomentum = pHelixFit->GetExtrapolatedMomentum(calorimeterPosition);
//     trackParameters.m_trackStateAtCalorimeter = pandora::TrackState(calorimeterPosition, calorimeterMomentum);

//     delete pHelixFit;
// }

