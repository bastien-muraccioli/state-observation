#ifndef HMMCONTACTESTIMATOR
#define HMMCONTACTESTIMATOR


#include"state-observation/observer/hmm.hpp"

#pragma once

#include<vector>
#include<string>

using namespace std;

namespace stateObservation{

    class STATE_OBSERVATION_DLLAPI hmm_contact_force : public hmm {

        public:

            hmm_contact_force(const int n, const int p, const std::vector<std::string>& labels, const hmm::Emission& emission, const hmm::Transition& transition);
            
            hmm_contact_force(int n, int p);

            void setMeasurement(Vector2 F, Vector3 v, Vector3 v_p, TimeIndex k);


    };

}
#endif
