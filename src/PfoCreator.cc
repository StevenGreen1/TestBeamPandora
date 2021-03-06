/**
 *  @file   MarlinPandora/src/PfoCreator.cc
 * 
 *  @brief  Implementation of the pfo creator class.
 * 
 *  $Log: $
 */

#include "marlin/Global.h"
#include "marlin/Processor.h"

#include "EVENT/LCCollection.h"

#include "IMPL/ClusterImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCFlagImpl.h"
#include "IMPL/ReconstructedParticleImpl.h"

#include "CalorimeterHitType.h"
#include "ClusterShapes.h"

#include "Api/PandoraApi.h"

#include "Objects/ParticleFlowObject.h"

#include "CalicePandoraProcessor.h"
#include "PfoCreator.h"

#include <cmath>

PfoCreator::PfoCreator(const Settings &settings) :
    m_settings(settings),
    m_pPandora(CalicePandoraProcessor::GetPandora())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

PfoCreator::~PfoCreator()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode PfoCreator::CreateParticleFlowObjects(EVENT::LCEvent *pLCEvent)
{
    pandora::ParticleFlowObjectList particleFlowObjectList;
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::GetParticleFlowObjects(*m_pPandora, particleFlowObjectList));

    IMPL::LCCollectionVec *pClusterCollection = new IMPL::LCCollectionVec(LCIO::CLUSTER);
    IMPL::LCCollectionVec *pReconstructedParticleCollection = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

    IMPL::LCFlagImpl lcFlagImpl(pClusterCollection->getFlag());
    lcFlagImpl.setBit(LCIO::CLBIT_HITS);
    pClusterCollection->setFlag(lcFlagImpl.getFlag());

    std::vector<std::string> subDetectorNames ;
    subDetectorNames.push_back("ecal") ; const unsigned int ecal_Index(0) ;
    subDetectorNames.push_back("hcal") ; const unsigned int hcal_Index(1) ;
    subDetectorNames.push_back("tc")   ; const unsigned int tc_Index(2)   ;

    pClusterCollection->parameters().setValues("ClusterSubdetectorNames", subDetectorNames);

    // Create lcio "reconstructed particles" from the pandora "particle flow objects"
    for (pandora::ParticleFlowObjectList::iterator itPFO = particleFlowObjectList.begin(), itPFOEnd = particleFlowObjectList.end();
         itPFO != itPFOEnd; ++itPFO)
    {
        IMPL::ReconstructedParticleImpl *pReconstructedParticle= new ReconstructedParticleImpl();

        pandora::ClusterAddressList clusterAddressList = (*itPFO)->GetClusterAddressList();
//         pandora::TrackAddressList trackAddressList = (*itPFO)->GetTrackAddressList();

        // Create lcio clusters
        for (pandora::ClusterAddressList::iterator itCluster = clusterAddressList.begin(), itClusterEnd = clusterAddressList.end();
            itCluster != itClusterEnd; ++itCluster)
        {
            IMPL::ClusterImpl *pCluster = new ClusterImpl();

            const unsigned int nHitsInCluster((*itCluster).size());

            float clusterEnergy(0.);
            float *pHitE = new float[nHitsInCluster];
            float *pHitX = new float[nHitsInCluster];
            float *pHitY = new float[nHitsInCluster];
            float *pHitZ = new float[nHitsInCluster];

            for (unsigned int iHit = 0; iHit < nHitsInCluster; ++iHit)
            {
                EVENT::CalorimeterHit *pCalorimeterHit = (CalorimeterHit*)((*itCluster)[iHit]);
                pCluster->addHit(pCalorimeterHit, 1.0);

                const float caloHitEnergy(pCalorimeterHit->getEnergy());
                clusterEnergy += caloHitEnergy;

                pHitE[iHit] = caloHitEnergy;
                pHitX[iHit] = pCalorimeterHit->getPosition()[0];
                pHitY[iHit] = pCalorimeterHit->getPosition()[1];
                pHitZ[iHit] = pCalorimeterHit->getPosition()[2];

//                 std::vector<float> &subDetectorEnergies = pCluster->subdetectorEnergies();
//                 subDetectorEnergies.resize(subDetectorNames.size());

//                 std::cout << "subdet type " << CHT(pCalorimeterHit->getType()).caloID() << std::endl;
//                 switch (CHT(pCalorimeterHit->getType()).caloID())
//                 {
//                 case CHT::ecal:  subDetectorEnergies[ecal_Index ] += caloHitEnergy; break;
//                 case CHT::hcal:  subDetectorEnergies[hcal_Index ] += caloHitEnergy; break;
//                 case CHT::yoke:  subDetectorEnergies[tc_Index   ] += caloHitEnergy; break;
//                 default: streamlog_out(DEBUG) << "no subdetector found for hit with type: " << pCalorimeterHit->getType() 
//                                               << "  ECAL is: " << CHT::ecal 
//                                               << "  HCAL is: " << CHT::hcal 
//                                               << "  tail catcher (yoke) is: " << CHT::yoke << std::endl;
//                 }
            }

            pCluster->setEnergy(clusterEnergy);

            ClusterShapes *pClusterShapes = new ClusterShapes(nHitsInCluster, pHitE, pHitX, pHitY, pHitZ);
            pCluster->setPosition(pClusterShapes->getCentreOfGravity());
            pCluster->setIPhi(std::atan2(pClusterShapes->getEigenVecInertia()[1], pClusterShapes->getEigenVecInertia()[0]));
            pCluster->setITheta(std::acos(pClusterShapes->getEigenVecInertia()[2]));

            pClusterCollection->addElement(pCluster);
            pReconstructedParticle->addCluster(pCluster);

            delete pClusterShapes;
            delete[] pHitE; delete[] pHitX; delete[] pHitY; delete[] pHitZ;
        }

//         // Add tracks to the lcio reconstructed particles
//         for (pandora::TrackAddressList::iterator itTrack = trackAddressList.begin(), itTrackEnd = trackAddressList.end(); itTrack != itTrackEnd;
//             ++itTrack)
//         {
//             pReconstructedParticle->addTrack((Track*)(*itTrack));
//         }

        float momentum[3] = {(*itPFO)->GetMomentum().GetX(), (*itPFO)->GetMomentum().GetY(), (*itPFO)->GetMomentum().GetZ()};
        pReconstructedParticle->setMomentum(momentum);
        pReconstructedParticle->setEnergy((*itPFO)->GetEnergy());
        pReconstructedParticle->setMass((*itPFO)->GetMass());
        pReconstructedParticle->setCharge((*itPFO)->GetCharge());
        pReconstructedParticle->setType((*itPFO)->GetParticleId());

        pReconstructedParticleCollection->addElement(pReconstructedParticle);
    }


    pLCEvent->addCollection(pClusterCollection, m_settings.m_clusterCollectionName.c_str());
    pLCEvent->addCollection(pReconstructedParticleCollection, m_settings.m_pfoCollectionName.c_str());

    return pandora::STATUS_CODE_SUCCESS;
}
