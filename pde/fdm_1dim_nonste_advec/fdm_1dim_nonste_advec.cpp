// 1dim non-steady advection equation by fdm
// use upwind difference method
// http://www.slis.tsukuba.ac.jp/~fujisawa.makoto.fu/lecture/mic/text/09_derivative2.pdf を参照

#include <cmath>    // for std::sin
#include <fstream>  // for std::ofstream 
#include <iostream> // for std::cerr, std::endl
#include <sstream>  // for std::string, std::stringstream
#include <vector>   // for std::vector
#include <boost/assert.hpp>                   // for BOOST_ASSERT
#include <boost/format.hpp>                   // for boost::format
#include <boost/math/constants/constants.hpp> // for boost::math::constants::pi

namespace{
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

    static auto constexpr XRANGE = 1.0;
    static auto constexpr MESH = 100;
    static auto constexpr NODE = MESH + 1;
    static auto constexpr DX = XRANGE / MESH; 

    // static auto constexpr COURANT = 0.25;
    // static auto constexpr DT = COURANT * DX / ADVEC;
    // static auto constexpr TEND = 1.0;
    // static auto constexpr TLOOP = TEND / DT;

    static auto constexpr TEND = 10.0;
    static auto constexpr TLOOP = 10000;
    static auto constexpr TREP = 100;
    static auto constexpr DT = TEND / TLOOP;

   static auto constexpr ADVEC = 1.0; // advection coefficient
    static auto constexpr COURANT = ADVEC * DT / DX; // Courant number

    bool result_output(int const t, myvec const &u0);

    myvec discrete();

    myvec initial(myvec const &x);

    void upwind_diffe(myvec &u0);
}

int main(){
    // COURANT > 1.0のとき強制終了 
    // CFL条件より, COURANT <= 1.0でないと正しい結果が得られない
    BOOST_ASSERT(COURANT <= 1.0);

    auto x = discrete();

    auto u0 = initial(x);

    if (!result_output(0, u0)){
		std::cerr << "output file not open" << std::endl;
        return -1;
    }

    for (auto t = 0; t <= TLOOP; t++){
        upwind_diffe(u0);

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

		for (auto i = 0; i < NODE; i++){
			std::fprintf(fp.get(), "%.2f %.7f\n", static_cast<double>(i) * DX, u0[i]);
		}

		return true;
    }

    void upwind_diffe(myvec &u0){
        myvec u(NODE, 0.0);

        switch (BCT){
            case boundary_conditon_type::DIRICLET:
                u[0] = D0;
                u[NODE - 1] = D1;
                break;

            case boundary_conditon_type::LEFT_NEUMANN:
                u[0] = u0[0] - COURANT * DX * N0;
                u[NODE - 1] = D1;
                break;

            case boundary_conditon_type::RIGHT_NEUMANN:
                u[0] = D0;
                u[NODE - 1] = u0[NODE - 1] - COURANT * DX * N1;
                break;

            case boundary_conditon_type::PERIODIC:
                u[0] = u0[0] - COURANT * 0.5 * (u0[1] - u0[NODE - 1]);
                u[NODE - 1] = u0[NODE - 1] - COURANT * 0.5 * (u0[0] - u0[NODE - 2]);
                break;

            default:
                BOOST_ASSERT(!"boundary_conditon_type is abnormal value");
                break;
        }

        for (auto i = 1; i < NODE - 1; i++){
            // upwind difference method
            u[i] = u0[i] - COURANT * 0.5 * (u0[i] - u0[i - 1]);
        }

        // update of u
        u0 = std::move(u);
    }

    myvec discrete(){
        myvec x(NODE, 0.0);

        // discrization of xrange
        for (auto i = 0; i < NODE; i++){
            x[i] = static_cast<double>(i) * DX;
        }

        return x;
    }

    myvec initial(myvec const &x){
        using namespace boost::math::constants; 
        myvec u0(NODE, 0.0);

        // initial condition of u
        for (auto i = 0; i < NODE; i++){
            u0[i] = 0.5 * (std::cos(2.0 * pi<double>() * x[i]) + 1);
        }

        return u0;
    }
}
