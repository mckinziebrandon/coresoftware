#include "RawClusterBuilderIA.h"

// User-defined includes (and boost).
#include "IslandAlgorithm.h" // Use this instead of PHMakeGroups.h
#include "IslandAlgorithmTower.h"

// G4CEMC includes.
#include "RawTower.h"
#include "RawTowerGeomContainer.h"
#include "RawTowerContainer.h"
#include "RawClusterv1.h"
#include "RawClusterContainer.h"

#define BOOST_NO_HASH // Our version of boost.graph is incompatible with GCC-4.3 w/o this flag
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
// Fun4All/PHENIX includes.
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

const std::string PATH = "~/bmckinz/RawClusterBuilderIA/rootFiles/";
typedef std::multimap<int, IslandAlgorithmTower>             TowerMap;
typedef std::pair<const int, IslandAlgorithmTower>           TowerPair;
typedef std::pair<const unsigned int, RawTower*> RawTowerPair;

/* ------------------------------------------------------ *
   RawClusterBuilderIA::RawClusterBuilderIA()
   1. Call parent constructor (SubsysReco).                     
   2. Initialize values of any desired attributes.
 * -------------------------------------------------------- */
RawClusterBuilderIA::RawClusterBuilderIA(const std::string& name)
    : SubsysReco(name),
      _clusters(NULL),
      _min_tower_e(0.0),
      _checkEnergyConserv(0),
      _clusterSimple(false), 
      _detector("NONE") {}

/* ----------------------------------------------------- *
   RawClusterBuilderIA::InitRun()
   1. Setup clusterNode and _clusters object. (CreateNodes)
   2. Create ROOT objects for file IO (temporary).
 * ----------------------------------------------------- */
int RawClusterBuilderIA::Init(PHCompositeNode* topNode) {
    try {
        _CreateNodes(topNode);
    } catch (std::exception &e) {
        std::cout << PHWHERE << ": " << e.what() << std::endl;
        throw;
    }

    std::cout << " - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - "     << std::endl;
    std::cout << "[RawClusterBuilderIA::Init] PARTICLE TYPE IS " << _particleType   << std::endl;
    std::cout << " - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - "   << std::endl;

    _fileName   = Form("%srcb_%s_%dGeV.root", PATH.c_str(), _particleType.c_str(), (int)(_genPT*10));
    _file       = new TFile(_fileName.c_str(),"RECREATE"); 
    gROOT->ProcessLine("#include <vector>");

    std::string varList;
    varList     = "iEvent:towerID:clusterID:energy:ET:eta:phi:ieta:iphi:nBinsEta:nBinsPhi";
    ntp_tower   = new TNtuple("ntp_tower", "tower values", varList.data());

    _tCluster = new TTree("tCluster", "cluster tree");
    _tCluster->Branch("towerIDs",  &_towerIDs);
    _tCluster->Branch("iEvent", &_iEvent);
    _tCluster->Branch("clusterID", &_clusterID);
    _tCluster->Branch("genPT", &_genPT);
    _tCluster->Branch("recoEnergy", &_energy);
    _tCluster->Branch("recoET", &_ET);
    _tCluster->Branch("eta", &_eta);
    _tCluster->Branch("phi", &_phi);
    _tCluster->Branch("nClusters", &_nClusters);
    _tCluster->Branch("nTowers", &_nTowers);

    return Fun4AllReturnCodes::EVENT_OK;
}

/* ------------------------------------------------------------ *
   RawClusterBuilderIA::process_event(...)
   1. Clear containers of any previous event info.
   2. Grab _towers and _towerGeom from the node tree.
   3. Obtain clusters from IslandAlgorithm functions.
   4. Store cluster information in _clusters object.
   5. Calculate cluster quantities (e.g. energy).
 * ------------------------------------------------------------ */
int RawClusterBuilderIA::process_event(PHCompositeNode *topNode) {
    namespace IAlgorithm = IslandAlgorithm;  
    _towerIDs.clear();
    _energyVec.clear();
    _ETVec.clear();
    _etaVec.clear();
    _phiVec.clear();
    std::string nodeName;

    // Grab the container of RawTowers (kinematic info).
    nodeName = "TOWER_CALIB_" + _detector;
    _towers  = findNode::getClass<RTContainer>(topNode, nodeName.c_str());
    if (!_towers) return _NodeError(nodeName, Fun4AllReturnCodes::DISCARDEVENT);

    // Grab the container of RawTowerGeoms (tower geometry info).
    nodeName    = "TOWERGEOM_" + _detector;
    _towerGeom  = findNode::getClass<RTGeomContainer>(topNode, nodeName.c_str());
    if (!_towerGeom) return _NodeError(nodeName, Fun4AllReturnCodes::ABORTEVENT);

    // Store the number of bins in phi as a static value; No need to repeatedly set.
    IslandAlgorithmTower::setMaxPhiBin(_towerGeom->get_phibins());
    IslandAlgorithmTower::setMaxEtaBin(_towerGeom->get_etabins());

    // ------------------------------------------------------------------------------------------
    // The Island Algorithm:

    set_threshold_energy(0.1);
    // 1. Construct list of seed towers, defined as having energy above some threshold. 
    std::list<IslandAlgorithmTower> seedTowers     = IAlgorithm::GetSeedTowers(_towers, _towerGeom, _min_tower_e);
    // 2. Cluster the towers via searching method along eta and phi from each seed.
    TowerMap clusteredTowers;
    if (_clusterSimple) clusteredTowers = IAlgorithm::GetSimpleClusters(seedTowers, _towers, _towerGeom);
    else                clusteredTowers = IAlgorithm::GetClusteredTowers(seedTowers, _towers, _towerGeom);

    // ------------------------------------------------------------------------------------------

    // Fill _clusters (now empty) with the clusteredTowers and calculate their values.
    foreach (TowerPair& towerPair, clusteredTowers) {
        // Either get the cluster with this ID or create a new one.
        RawCluster* rawCluster = _clusters->getCluster(towerPair.first); 
        if (!rawCluster) _CreateNewCluster(&rawCluster);
        // Add the tower to this cluster.
        rawCluster->addTower(towerPair.second.getID(), towerPair.second.getEnergy());
        if (verbosity) _PrintCluster(towerPair);
    }

    std::cout << "seedTowers.size() = "      << seedTowers.size()        << std::endl;
    std::cout << "clusteredTowers.size() = " << clusteredTowers.size()   << std::endl;
    std::cout << "_clusters->size() = "      << _clusters->size()        << std::endl;

    // Calculate/store energy, eta, phi of clusters given clusteredTowers information.
    _FillClustersEnergy(clusteredTowers);
    _FillClustersEta(clusteredTowers);
    _FillClustersPhi(clusteredTowers);
    for (unsigned i = 0; i < _clusters->size(); i++) {
        _AssignClusterValues(i);
    }

    // Handle ROOT I/O. 
    //_FillTowerTree(_GetAllTowers());
    _FillTowerTree(clusteredTowers);
    _FillClusterTree();

    if (_checkEnergyConserv) _CheckEnergyConservation();
    return Fun4AllReturnCodes::EVENT_OK;

}

int RawClusterBuilderIA::End(PHCompositeNode *topNode) {
    _ShowTreeEntries();
    std::cout << "WRITING TO FILE " << _file->GetName() << std::endl;
    _file->Write();
    _file->Close();
    return Fun4AllReturnCodes::EVENT_OK;
}

// -----------------------------------------------------------------------------
// -------------------------- PRIVATE HELPER METHODS. --------------------------
// -----------------------------------------------------------------------------

void RawClusterBuilderIA::_AssignClusterValues(int iCluster) {
    RawCluster* cluster = _clusters->getCluster(iCluster);
    cluster->set_energy(_ETVec[iCluster]);
    cluster->set_eta(_etaVec[iCluster]);
    cluster->set_phi(_phiVec[iCluster]);
    if (verbosity) {
        std::cout << " (eta,phi,e) = (" << cluster->get_eta() << ", "
            << cluster->get_phi() << ","
            << cluster->get_energy() << ")"
            << std::endl;
    }
}

// Return list of all RawTower pairs in _towers->getTowers() converted to IslandAlgorithmTowers.
std::list<IslandAlgorithmTower> RawClusterBuilderIA::_GetAllTowers() {
    std::list<IslandAlgorithmTower> allTowers;
    foreach (RawTowerPair& towerPair, _towers->getTowers()) {
        // TODO : change order of arguments. 
        _InsertTower(allTowers, towerPair);
    }
    return allTowers;
}
// Given iterator to a seed tower, place relevant info into std::vector of seed towers.
void RawClusterBuilderIA::_InsertTower(std::list<IslandAlgorithmTower>&  towerList, RawTowerPair towerPair)  {
    IslandAlgorithmTower rtHelper(towerPair.second);
    rtHelper.setCenter(_towerGeom);
    towerList.push_back(rtHelper);
}

void RawClusterBuilderIA::_FillClustersEnergy(TowerMap clusteredTowers) {
    foreach (TowerPair& towerPair, clusteredTowers) {
        IslandAlgorithmTower tower = towerPair.second;
        _energyVec[towerPair.first] += tower.getEnergy();  
        _ETVec[towerPair.first]     += tower.getET();
    }
}


void RawClusterBuilderIA::_FillClustersEta(TowerMap clusteredTowers) {
    foreach (TowerPair& towerPair, clusteredTowers) {
        IslandAlgorithmTower tower = towerPair.second;
        _etaVec[towerPair.first] += tower.getET() * tower.getEtaCenter();
    }
    for (unsigned i = 0; i < _clusters->size(); i++) {
        _etaVec[i] = (_ETVec[i] > 0) ? _etaVec[i] / _ETVec[i] : 0.0;
    }
}

void RawClusterBuilderIA::_FillClustersPhi(TowerMap clusteredTowers) {
    foreach (TowerPair& towerPair, clusteredTowers) {
        IslandAlgorithmTower tower = towerPair.second;
        _phiVec[towerPair.first] += tower.getET() * tower.getPhiCenter();
    }
    // First, get all constitutent tower phi's as an energy-weighted sum.
    // Then divide by the total cluster energy.
    for (unsigned int i = 0; i < _clusters->size(); i++) {
        _phiVec[i] = (_ETVec[i] > 0) ? _phiVec[i] / _ETVec[i] : 0.0;
        if (_phiVec[i] > M_PI)  _phiVec[i] -= 2. * M_PI;
    }
    // Finally, correct the mean Phi calculation for clusters at Phi discontinuity.
    for (unsigned int iCluster = 0; iCluster < _clusters->size(); iCluster++) {
        RawCluster *cluster = _clusters->getCluster(iCluster);
        float oldPhi        = cluster->get_phi();
        bool corr           = _CorrectPhi(cluster);
        if (corr) {
            std::cout << PHWHERE << " Cluster Phi corrected: " 
                << oldPhi  << " " << cluster->get_phi() << std::endl;
        }
    }
}

// ----------------------------------------------------- *
// RawClusterBuilderIA::CreateNodes()                          *
// ----------------------------------------------------- //
void RawClusterBuilderIA::_CreateNodes(PHCompositeNode *topNode) {
    PHNodeIterator iter(topNode);
    // Grab the cEMC node
    PHCompositeNode *dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode) {
        std::cerr << PHWHERE << "DST Node missing, doing nothing." << std::endl;
        throw std::runtime_error("Failed to find DST node in EmcRawTowerBuilder::CreateNodes");
    }
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode*>(dstiter.findFirst("PHCompositeNode",_detector ));
    if(!DetNode){
        DetNode = new PHCompositeNode(_detector);
        dstNode->addNode(DetNode);
    }

    _clusters       = new RawClusterContainer();
    std::string nodeName = "CLUSTER_" + _detector;
    PHIODataNode<PHObject> *clusterNode = new PHIODataNode<PHObject>(_clusters, nodeName.c_str(), "PHObject");
    DetNode->addNode(clusterNode);
}


bool RawClusterBuilderIA::_CorrectPhi(RawCluster* cluster) {
    double sum      = cluster->get_energy();
    double phimin   = 999.;
    double phimax   = -999.;
    RawCluster::TowerConstRange begin_end = cluster->get_towers();
    RawCluster::TowerConstIterator iter;
    for (iter = begin_end.first; iter != begin_end.second; ++iter) { 
        RawTower* tmpt = _towers->getTower(iter->first);
        double phi = _towerGeom->get_phicenter(tmpt->get_binphi());
        if(phi > M_PI) phi = phi - 2.*M_PI; 
        if (phi < phimin) {
            phimin = phi;
        }
        if (phi > phimax) {
            phimax = phi;
        }
    }

    // cluster is not at phi discontinuity
    if ((phimax - phimin) < 3.) return false; 

    float mean = 0.;
    for (iter =begin_end.first; iter != begin_end.second; ++iter) { 
        RawTower* tmpt = _towers->getTower(iter->first);
        double e = tmpt->get_energy();
        double phi = _towerGeom->get_phicenter(tmpt->get_binphi());
        if(phi > M_PI) phi = phi - 2.*M_PI; 
        if (phi < 0.) {
            phi = phi + 2.*M_PI;  // shift phi range for correct mean calculation
        }
        mean += e * phi;
    }
    mean = mean / sum;
    if (mean > M_PI) {
        mean = mean - 2.*M_PI;  // shift back
    }
    cluster->set_phi(mean);
    return true; // mean phi was corrected
}
void RawClusterBuilderIA::_CheckEnergyConservation() {
    double ecluster = _clusters->getTotalEdep();
    double etower   = _towers->getTotalEdep();
    if (ecluster > 0 && (fabs(etower - ecluster) / ecluster) > 1e-9) {
        std::cout << "energy conservation violation: ETower: " << etower
            << " ECluster: " << ecluster 
            << " diff: " << etower - ecluster << std::endl;
    } else if (etower != 0) {
        std::cout << "energy conservation violation: ETower: " << etower
            << " ECluster: " << ecluster << std::endl;
    }
}

void RawClusterBuilderIA::_CreateNewCluster(RawCluster** rawCluster) {
    *rawCluster = new RawClusterv1();
    _clusters->AddCluster(*rawCluster);
    _energyVec.push_back(0.0);
    _ETVec.push_back(0.0);
    _etaVec.push_back(0.0);
    _phiVec.push_back(0.0);
}

// Serves to make ugly code in process_event less ugly.
int RawClusterBuilderIA::_NodeError(std::string nodeName, int retCode) {
    std::cout << PHWHERE << ": Could not find node " 
        << nodeName.data() << std::endl;
    return retCode;
}

void RawClusterBuilderIA::_PrintCluster(TowerPair towerPair) {
    std::cout << "RawClusterBuilderIA id: " << (towerPair.first)
        << " Tower: " << " (iEta,iPhi) = (" << towerPair.second.getBinEta() 
        << "," << towerPair.second.getBinPhi() << ") " << " (eta,phi,e) = (" 
        << towerPair.second.getEtaCenter() << ","
        << towerPair.second.getPhiCenter() << ","
        << towerPair.second.getEnergy() << ")"
        << std::endl;
}

// Functions related to ROOT I/O (temporary). 
void RawClusterBuilderIA::_ShowTreeEntries() {
    for (int i = 0; i < _tCluster->GetEntries(); i++) {
        _tCluster->Show(i);
    } 
}

void RawClusterBuilderIA::_FillTowerTree(TowerMap clusteredTowers) {
    foreach (TowerPair& towerPair, clusteredTowers) {
        int clusterID = towerPair.first;
        IslandAlgorithmTower tower = towerPair.second;
        ntp_tower->Fill(
                _iEvent,
                tower.getID(), 
                clusterID, 
                tower.getEnergy(), 
                tower.getET(),
                tower.getEtaCenter(), 
                tower.getPhiCenter(), 
                tower.getBinEta(), 
                tower.getBinPhi(), 
                tower.getMaxEtaBin(), 
                tower.getMaxPhiBin());
    }
}


void RawClusterBuilderIA::_FillTowerTree(std::list<IslandAlgorithmTower> allTowers) {
    foreach (IslandAlgorithmTower& tower, allTowers) {
        ntp_tower->Fill(
                tower.getID(), 
                tower.getEnergy(), 
                tower.getET(),
                tower.getEtaCenter(), 
                tower.getPhiCenter(), 
                tower.getBinEta(), 
                tower.getBinPhi(), 
                tower.getMaxEtaBin(), 
                tower.getMaxPhiBin(), 
                allTowers.size(), 
                _iEvent);

    }
}

void RawClusterBuilderIA::_FillClusterTree() {
    typedef std::pair<RawTowerDefs::keytype, float> TowIDEnergy;
    for (unsigned i = 0; i < _clusters->size(); i++) {

        RawCluster* rawCluster = _clusters->getCluster(i);

        foreach (TowIDEnergy towIDEnergy, rawCluster->get_towers()) {
            _towerIDs.push_back((int) towIDEnergy.first);
        }

        _clusterID = rawCluster->get_id();
        _energy = _energyVec[i];
        _ET = _ETVec[i];
        _eta = _etaVec[i];
        _phi = _phiVec[i];
        _nClusters = _clusters->size();
        _nTowers = rawCluster->getNTowers();

        _tCluster->Fill();
    }
}

