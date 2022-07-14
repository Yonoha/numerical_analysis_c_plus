// 2dim Lotka Volterra equation by boost::odeint

#include <array>    // for std::array
#include <fstream>  // for std::ofstream
#include <iostream> // for std::cerr, std::endl
#include <boost/format.hpp>         // for boost::format
#include <boost/numeric/odeint.hpp> // for boost::numeric::odeint

namespace{
    using state = std::array<double, 2>;

    static auto constexpr DT = 0.05;
    static auto constexpr T0 = 0.0;
    static auto constexpr TMAX = 5.0;

    class Lotka_Volterra_system final{
        double const alpha_;
        double const beta_;
        double const delta_;
        double const gamma_;

        public:
            Lotka_Volterra_system(double alpha, double beta, double delta, double gamma)
                : alpha_(alpha), beta_(beta), delta_(delta), gamma_(gamma)
            {}
            Lotka_Volterra_system() = delete;
            ~Lotka_Volterra_system() = default;

            void operator()(state const &x, state &dx, double t) const; //関数内でalpha等を変更しないためのconst
    }; 
}

int main(){   
    Lotka_Volterra_system System(2.0, 3.0, 5.0, 4.0);
    state State = {1.0, 0.5};

    //6th-Adams Bashforth Moulton method
    boost::numeric::odeint::adams_bashforth_moulton<6, state> Stepper;

    // open output file
    std::ofstream ofs("data_Lotka_Volterra.txt");

    if (!ofs){
        std::cerr << "output file not open" << std::endl;
        return -1;
    }

    //time = T0からtime = TMAXまで時間発展をDT刻みで計算
    boost::numeric::odeint::integrate_const(
        Stepper,
        System,
        State,
        T0,
        TMAX,
        DT,
        [&ofs](state const &x, double t) {ofs << boost::format("%.2f %.7f %.7f\n") % t % x[0] % x[1];}
    );
}

namespace{
    void Lotka_Volterra_system::operator()(state const &x, state &dx, double t) const{
        dx[0] = x[0] * (alpha_ - beta_ * x[1]);
        dx[1] = - x[1] * (gamma_ - delta_ * x[0]);
    }
}
