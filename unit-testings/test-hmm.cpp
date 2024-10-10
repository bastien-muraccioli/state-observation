#include<iostream>
#include<state-observation/observer/hmm-contact.hpp>

#include<vector>
#include<string>
#include<map>
#include<functional>

int main () {

    std::cout<<"Starting HMM test"<<std::endl;

    ////////////// HMM Contact ////////////

    std::cout<< "Test HMM Contact estimator" << std::endl;
       
    int n = 3;
    int p = 3;
    stateObservation::hmm_contact_force HMM(n, p);
    std::cout << "Test contruction passed !" << std::endl;

    stateObservation::Vector2 F = {1, 2};
    stateObservation::Vector3 v = {1, 2, 3};
    stateObservation::Vector3 v_p = {1, 2, 3};

    HMM.setMeasurement(F, v, v_p, 3);
    std::cout << "SetMeasurment test passed !" << std::endl;

    
    std::vector<std::string> labels = {"A", "B", "C"};
    stateObservation::hmm_contact_force::Emission emission;
    emission["A"] = [](stateObservation::hmm_contact_force::MeasureVector O_t){return 2*O_t(1);};
    emission["B"] = [](stateObservation::hmm_contact_force::MeasureVector O_t){return 2*O_t(1);};
    emission["C"] = [](stateObservation::hmm_contact_force::MeasureVector O_t){return 2*O_t(1);};

    stateObservation::hmm_contact_force::Transition transition;
    transition["A"]["B"] = [](stateObservation::hmm_contact_force::MeasureVector O_t){return 2*O_t(1);};
    transition["A"]["C"] = [](stateObservation::hmm_contact_force::MeasureVector O_t){return 2*O_t(1);};
    transition["A"]["A"] = [&transition](stateObservation::hmm_contact_force::MeasureVector O_t){return 1-transition["A"]["B"](O_t)-transition["A"]["C"](O_t);};

    transition["B"]["C"] = [](stateObservation::hmm_contact_force::MeasureVector O_t){return 2*O_t(1);};
    transition["B"]["A"] = [](stateObservation::hmm_contact_force::MeasureVector O_t){return 2*O_t(1);};
    transition["B"]["B"] = [&transition](stateObservation::hmm_contact_force::MeasureVector O_t){return 1-transition["B"]["C"](O_t)-transition["B"]["A"](O_t);};

    transition["C"]["A"] = [](stateObservation::hmm_contact_force::MeasureVector O_t){return 2*O_t(1);};
    transition["C"]["B"] = [](stateObservation::hmm_contact_force::MeasureVector O_t){return 2*O_t(1);};
    transition["C"]["C"] = [&transition](stateObservation::hmm_contact_force::MeasureVector O_t){return 1-transition["C"]["A"](O_t)-transition["C"]["B"](O_t);};


    HMM.setEmTr(labels, emission, transition);


    std::cout << "SetEmTr test passed !" << std::endl;
    return 0;
}
