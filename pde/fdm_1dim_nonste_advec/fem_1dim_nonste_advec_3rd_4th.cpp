// 1dim non-steady advection equation by fdm
// use 3rd accuracy upwind difference method
// use 4th accuracy Runge Kutta method

// http://www.slis.tsukuba.ac.jp/~fujisawa.makoto.fu/lecture/mic/text/09_derivative2.pdf を参照
// http://pengubeforenitis.g1.xrea.com/study/note/upwind_difference.pdf を参照

#include <cmath>    // for initial condition of u (if use)
#include <fstream>  // for std::ofstream 
#include <iostream> // for std::cerr, std::endl
#include <tuple>    // for std::tuple
#include <vector>   // for std::vector
#include <boost/assert.hpp>                   // for BOOST_ASSERT
#include <boost/format.hpp>                   // for boost::format
#include <boost/math/constants/constants.hpp> // for boost::math::constants::pi

namespace{
    using mytuple = std::tuple<int, int, int, int>;
    using myvec = std::vector<double>;

    // boundary condition
    enum class boundary_conditon_type{
        DIRICLET = 0,
        LEFT_NEUMANN = 1,
        RIGHT_NEUMANN = 2,
        PERIODIC = 3,
    };

    static auto constexpr BCT = boundary_conditon_type::PERIODIC;
    static auto constexpr D0 = 1.0; // left Diriclet
    static auto constexpr D1 = 1.0; // right Diriclet
    static auto constexpr N0 = 0.0; // left Neumann
    static auto constexpr N1 = 0.0; // right Neumann

    static auto constexpr TEND = 1.0;
    static auto constexpr TLOOP = 10000;
    static auto constexpr TREP = 100;
    static auto constexpr DT = TEND / TLOOP;

    static auto constexpr XRANGE = 1.0;
    static auto constexpr XMESH = 100;
    static auto constexpr XNODE = XMESH + 1;
    static auto constexpr DX = XRANGE / XMESH; 

    static auto constexpr ADVEC = -1.0; // advection coefficient
    static auto constexpr COURANT = ADVEC * DT / DX; // Courant number

    bool result_output(int const t, myvec const &u0);

    mytuple make_param(int i);

    myvec Runge_Kutta(myvec const u0, myvec const k);

    myvec upwind_diffe(myvec const ubefore);

    myvec uinitial(myvec const &x);

    myvec xdiscrete();

    void diffe_method(myvec &u0);
}

int main(){
    // std::fabs(COURANT) > 1.0のとき強制終了 
    // CFL条件より, std::fabs(COURANT) <= 1.0でないと正しい結果が得られない
    BOOST_ASSERT(std::fabs(COURANT) <= 1.0);

    auto x = xdiscrete();

    auto u0 = uinitial(x);

    if (!result_output(0, u0)){
		std::cerr << "output file not open" << std::endl;
        return -1;
    }

    for (auto t = 1; t <= TLOOP; t++){
        diffe_method(u0);

        if(!(t % (TLOOP / TREP))){
            if (!result_output(t, u0)){
		        std::cerr << "output file not open" << std::endl;
                return -1;
            }
        }
    }

    return 0;
}

namespace{
    bool result_output(int const t, myvec const &u0){
		auto const filename = boost::format("data_fdm_1dim_nonste_advec_%d.txt") % (static_cast<int>(t / (TLOOP / TREP)));

		std::unique_ptr<FILE, decltype(&std::fclose)> fp(std::fopen(filename.str().c_str(), "w"), std::fclose);

		if (!fp){
			return false;
		}

		for (auto i = 0; i < XNODE; i++){
			std::fprintf(fp.get(), "%.2f %.7f\n", static_cast<double>(i) * DX, u0[i]);
		}

		return true;
    }

    void diffe_method(myvec &u0){
        myvec k1(XNODE, 0.0);
        myvec k2(XNODE, 0.0);
        myvec k3(XNODE, 0.0);
        myvec k4(XNODE, 0.0);

        myvec u1(XNODE, 0.0);
        myvec u2(XNODE, 0.0);
        myvec u3(XNODE, 0.0);
        myvec u(XNODE, 0.0);

        switch (BCT){
            case boundary_conditon_type::DIRICLET:
                break;

            case boundary_conditon_type::LEFT_NEUMANN:
                break;

            case boundary_conditon_type::RIGHT_NEUMANN:
                break;

            case boundary_conditon_type::PERIODIC:
                for (auto i = 0; i < XNODE; i++){
                    k1 = upwind_diffe(u0);
                    u1 = Runge_Kutta(u0, k1);
                }

                for (auto i = 0; i < XNODE; i++){
                    k2 = upwind_diffe(u1);
                    u2 = Runge_Kutta(u0, k2);
                }

                for (auto i = 0; i < XNODE; i++){
                    k3 = upwind_diffe(u2);
                    u3 = Runge_Kutta(u0, k3);
                }

                for (auto i = 0; i < XNODE; i++){
                    k4 = upwind_diffe(u3);
                    u[i] = u0[i] - (COURANT / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) * DX;
                }
                break;

            default:
                BOOST_ASSERT(!"boundary_conditon_type is abnormal value");
                break;
        }

        // update of u
        u0 = std::move(u);
    }

    mytuple make_param(int i){
        auto left1 = i - 1;
        auto left2 = i - 2;
        auto right1 = i + 1;
        auto right2 = i + 2;

        if (i == 0){
            left1 = XMESH - 1;
            left2 = XMESH - 2;
        }
        else if (i == 1){
            left2 = XMESH - 1;
        }

        else if (i == XMESH - 1){
            right2 = 0;
        }

        else if (i == XMESH){
            right1 = 0;
            right2 = 1;
        }

        return std::make_tuple(left1, left2, right1, right2);
    }

    myvec Runge_Kutta(myvec const u0, myvec const k){
        myvec uafter(XNODE, 0.0);

        for (auto i = 0; i < XNODE; i++){
            uafter[i] = u0[i] - COURANT * k[i] * DX / 2.0;
        }

        return uafter;
    }
    
    myvec upwind_diffe(myvec const ubefore){
        myvec k(XNODE, 0.0);

        for (auto i = 0; i < XNODE; i++){
            auto const [left1, left2, right1, right2] = make_param(i);

            if (ADVEC >= 0.0){
                k[i] = (2.0 * ubefore[right1] + 3.0 * ubefore[i] - 6.0 * ubefore[left1] + ubefore[left2]) / (6.0 * DX);
            }

            else{
                k[i] = (-2.0 * ubefore[left1] - 3.0 * ubefore[i] + 6.0 * ubefore[right1] - ubefore[right2]) / (6.0 * DX);
            }
        } 

        return k;     
    }
    
    myvec uinitial(myvec const &x){
        using namespace boost::math::constants; 

        myvec u0(XNODE, 0.0);

        // initial condition of u
        for (auto i = 0; i < XNODE; i++){
            u0[i] = 0.5 * (std::cos(2.0 * pi<double>() * x[i]) + 1);
        }

        return u0;
    }

    myvec xdiscrete(){
        myvec x(XNODE, 0.0);

        // discrization of xrange
        for (auto i = 0; i < XNODE; i++){
            x[i] = static_cast<double>(i) * DX;
        }

        return x;
    }
}
