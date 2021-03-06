#ifndef __TRUTHJETINPUT_H__
#define __TRUTHJETINPUT_H__

#include <phool/PHCompositeNode.h>

#include "JetInput.h"
#include "Jet.h"

#include <vector>

class TruthJetInput : public JetInput {
  
public:

  TruthJetInput(Jet::SRC input);
  virtual ~TruthJetInput() {}

  //! by default, TruthJetInput process all truth primary particle.
  //! However, it can be configured to read only one or more embedded stream via add_embedding_flag()
  //! It can be useful for reconstruct truth jet for embedded pythia jets only, etc.
  //! Call add_embedding_flag() multiple times to add multiple embed stream
  void add_embedding_flag(const int embed_stream_id)
  {_embed_id.push_back(embed_stream_id);}

  void identify(std::ostream& os = std::cout);
  
  Jet::SRC get_src() {return _input;}
  
  std::vector<Jet*> get_input(PHCompositeNode *topNode);
  
  void set_eta_range(float eta_min, float eta_max) {
    _eta_min = eta_min;
    _eta_max = eta_max;
  }
    
private:
  int _verbosity;
  Jet::SRC _input;
  float _eta_min;
  float _eta_max;

  //! if empty: process all primary particles
  //! if non-empty: only process primary particles in the selected embed stream.
  std::vector<int> _embed_id;

  bool use_embed_stream() {return _embed_id.size()>0;}
};

#endif
