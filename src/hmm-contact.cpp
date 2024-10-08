using namespace std;

#include "state-observation/observer/hmm-contact.hpp"


#include<string>
#include<map>
#include<functional>

namespace stateObservation {





    hmm_contact_force::hmm_contact_force(const int n, const int m, const int p, const std::vector<std::string>& labels, const Emission& emission, const Transition& transition): hmm{n, m, p, labels, emission, transition}{}
    
    hmm_contact_force::hmm_contact_force(int n, int m, int p): hmm{n, m, p}{}


    void hmm_contact_force::setMeasurement(Vector2 F, Vector3 v, Vector3 v_p, TimeIndex k){
        ObserverBase::MeasureVector y(8);
        y << F,v,v_p;
        ZeroDelayObserver::setMeasurement(y, k);
    }

    void hmm_contact_force::init(vector<string> labels, Emission emission, Transition transition){
        for(size_t i=0; i<hmm::n_state; i++){
            hmm::setEmTr(labels, emission, transition);
        }
    }

}

