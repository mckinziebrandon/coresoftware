#ifndef __RawClusterBuilderIA_H__
#define __RawClusterBuilderIA_H__

// C++ includes.
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>
#include <list>
// Fun4All/PHENIX includes.
#include <fun4all/SubsysReco.h>
// ROOT Includes.
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TROOT.h"

// Forward class declarations.
class RawTower;
class PHCompositeNode;
class RawCluster;
class RawClusterContainer;
class RawTowerContainer;
class RawTowerGeomContainer;
class IslandAlgorithmTower;

typedef RawTowerContainer       RTContainer;
typedef RawTowerGeomContainer   RTGeomContainer;
// Commenting these out before pull request. Appears that sPHENIX has similarly named typedefs. 
//typedef std::multimap<int, IslandAlgorithmTower>             TowerMap;
//typedef std::pair<const int, IslandAlgorithmTower>           TowerPair;
//typedef std::pair<const unsigned int, RawTower*> RawTowerPair;

class RawClusterBuilderIA : public SubsysReco {
    public:
        RawClusterBuilderIA(const std::string& name = "RawClusterBuilderIA");
        virtual ~RawClusterBuilderIA() {}
        int Init(PHCompositeNode *topNode);
        int process_event(PHCompositeNode *topNode);
        int End(PHCompositeNode *topNode);
        void Detector(const std::string &d)         { _detector = d; }
        void set_threshold_energy(const float e)    { _min_tower_e = e; }
        void checkenergy(const int i = 1)           { _checkEnergyConserv = i; }
        void SetGenPT(float pt)                     { _genPT = pt; }
        void SetParticleType(std::string s)         { _particleType = s; }
        void SetEvent(int i)                        { _iEvent = i; }
        void ClusterSimple(bool b)                  { _clusterSimple = b; }
    private:
        // Variables in constructor initializer list.
        RawClusterContainer* _clusters;     // default = NULL
        float   _min_tower_e;               // default = 0.0
        int     _checkEnergyConserv;        // default = 0
        bool    _clusterSimple;             // default = false
        std::string _detector;              // default = "NONE"

        // Other sPHENIX-type private variables.
        RTContainer*        _towers;
        RTGeomContainer*    _towerGeom;

        // Containers for cluster vars. 
        std::vector<float>  _energyVec; 
        std::vector<float>  _ETVec; 
        std::vector<float>  _etaVec; 
        std::vector<float>  _phiVec;

        // Private helper methods. 
        void _AssignClusterValues(int iCluster);
        void _CheckEnergyConservation();
        bool _CorrectPhi(RawCluster*);
        void _CreateNewCluster(RawCluster**);
        void _CreateNodes(PHCompositeNode *topNode);
        void _FillClustersEnergy(std::multimap<int, IslandAlgorithmTower>);
        void _FillClustersEta(std::multimap<int, IslandAlgorithmTower>);
        void _FillClustersPhi(std::multimap<int, IslandAlgorithmTower>);
        std::list<IslandAlgorithmTower> _GetAllTowers();
        void _InsertTower(std::list<IslandAlgorithmTower>&, std::pair<const unsigned, RawTower*>);
        int  _NodeError(std::string nodeName, int retCode);
        void _PrintCluster(std::pair<const int, IslandAlgorithmTower>);

        // _____ ROOT I/O Objects. (Temporary) _____
        std::string  _particleType; // for distinguishing single-events from particle gun.
        std::string _fileName;
        TFile*   _file;
        TNtuple* ntp_tower;

        // Cluster tree with (in order) branch variables.
        TTree* _tCluster;
        std::vector<int> _towerIDs;
        int     _iEvent;
        int     _clusterID;
        float   _genPT;
        float   _energy;
        float   _ET;
        float   _eta;
        float   _phi;
        int     _nClusters;
        int     _nTowers;

        void _FillTowerTree(std::list<IslandAlgorithmTower>);
        void _FillTowerTree(std::multimap<int,IslandAlgorithmTower>);
        void _FillClusterTree();
        void _ShowTreeEntries();

};

#endif
