// 1dim non-steady diffusion equation by fdm
// use center differential method

#include <cmath>    // for std::sin
#include <cstdio>   // for std::fclose, std::fopen, std::fprintf
#include <iostream> // for std::endl
#include <memory>   // for std::unique_ptr
#include <vector>   // fos std::vector
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

    static auto constexpr D0 = 0.0; // left Diriclet
    static auto constexpr D1 = 0.0; // right Diriclet
    static auto constexpr N0 = 0.0; // left Neumann
    static auto constexpr N1 = 0.0; // right Neumann

    static auto constexpr TEND = 1.0;
    static auto constexpr TLOOP = 10000;
    static auto constexpr TREP = 100;
    static auto constexpr DT = TEND / TLOOP;

    static auto constexpr XRANGE = 1.0;
    static auto constexpr MESH = 100;
    static auto constexpr NODE = MESH + 1;
    static auto constexpr DX = XRANGE / MESH; 

    // coefficient
    static auto constexpr DIFF = 0.01; // diffusion number
    static auto constexpr KAP = DIFF * DT / (DX * DX);

    bool result_output(int const t, myvec const &u0);

    myvec discrete();

    myvec initial(myvec const &x);

    void center_diffe(myvec &u0);
}

int main(){
    // KAP > 0.5のとき強制終了 
    // Neumannの安定性解析より, KAP <= 0.5でないと正しい結果が得られない
    BOOST_ASSERT(KAP <= 0.5); 

    auto x = discrete();

    auto u0 = initial(x);

    if (!result_output(0, u0)){
		std::cerr << "output file not open" << std::endl;
        return -1;
    }

    for (auto t = 0; t <= TLOOP; t++){
        center_diffe(u0);

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
		auto const filename = boost::format("data_fdm_1dim_nonste_diff_%d.txt") % (static_cast<int>(t / (TLOOP / TREP)));

		std::unique_ptr<FILE, decltype(&std::fclose)> fp(std::fopen(filename.str().c_str(), "w"), std::fclose);

		if (!fp){
			return false;
		}

		for (auto i = 0; i < NODE; i++){
			std::fprintf(fp.get(), "%.2f %.7f\n", static_cast<double>(i) * DX, u0[i]);
		}

		return true;
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

        // initial condition 
        for (auto i = 0; i < NODE; i++){
            u0[i] = std::sin(pi<double>() * x[i]);
        }

        return u0;
    }

    void center_diffe(myvec &u0){
        myvec u(NODE, 0.0);

        switch (BCT){
            case boundary_conditon_type::DIRICLET:
                u[0] = D0;
                u[NODE - 1] = D1;
                break;

            case boundary_conditon_type::LEFT_NEUMANN:
                u[0] = u0[0] + KAP * (u0[1] - u0[0] - DX * N0);
                u[NODE - 1] = D1;
                break;

            case boundary_conditon_type::RIGHT_NEUMANN:
                u[0] = D0;
                u[NODE - 1] = u0[NODE - 1] + KAP * (DX * N1 - u0[NODE - 1] + u0[NODE - 2]);
                break;

            case boundary_conditon_type::PERIODIC:
                u[0] = u0[0] + KAP * (u0[1] - 2.0 * u0[0] + u0[NODE - 1]);
                // u[NODE - 1] = u[0];
                u[NODE - 1] = u0[0] + KAP * (u0[0] - 2.0 * u0[NODE - 1] + u0[NODE - 2]);
                break;
        }

        for (auto j = 1; j < NODE - 1; j++){
            u[j] = u0[j] + KAP * (u0[j + 1] - 2.0 * u0[j] + u0[j - 1]);
        }

        // update of u
        for (auto j = 0; j < NODE; j++){
            u0[j] = u[j];
        }
    }
}
