// 2dim non-steady diffusion equation by fdm

#include <cmath>    // for initial condition (if use)
#include <cstdio>   // for std::fclose, std::fopen, std::fprintf
#include <iostream> // for std::cerr, std::endl
#include <memory>   // for std::unique_ptr
#include <utility>  // for std::move
#include <vector>   // for std::vector
#include <boost/assert.hpp>                   // for BOOST_ASSERT
#include <boost/format.hpp>                   // for boost::format
#include <boost/math/constants/constants.hpp> // for boost::math::constants::pi
#include <boost/multi_array.hpp>              // for boost::multi_array

namespace{
    using myvec = std::vector<double>;
    using my2darray = boost::multi_array<double, 2>;

    // boundary condition
    enum class boundary_conditon_type{
        DIRICLET = 0,
        LEFT_NEUMANN = 1,
        RIGHT_NEUMANN = 2,
        PERIODIC = 3,
    };
    static auto constexpr BCT = boundary_conditon_type::DIRICLET;

    static auto constexpr TEND = 10.0;
    static auto constexpr TREP = 100;
    static auto constexpr TLOOP = 10000; // TREP loopごとにgifにするのでTLOOP > 10 * TREP程度が良い
    static auto constexpr DT = TEND / TLOOP;

    static auto constexpr XRANGE = 2.0;
    static auto constexpr XMESH = 200;
    static auto constexpr XNODE = XMESH + 1;
    static auto constexpr DX = XRANGE / XMESH; 

    static auto constexpr YRANGE = 3.0;
    static auto constexpr YMESH = 300;
    static auto constexpr YNODE = YMESH + 1;
    static auto constexpr DY = YRANGE / YMESH;

    static auto constexpr DIFF = 0.01; // diffusion number
    static auto constexpr XKAP = DIFF * DT / (DX * DX);
    static auto constexpr YKAP = DIFF * DT / (DY * DY);
    static auto constexpr KAP = XKAP + YKAP;

    // output result for each i
    bool result_output(int const i, my2darray const &u0);

    // discrezation of x and y
    myvec xdiscrete();
    myvec ydiscrete();

    // initial condition of u(x, y, t)
    my2darray uinitial(myvec const &x, myvec const &y);

    // center differential method
    void center_diffe(my2darray &u0);
}

int main(){
    // KAP > 0.5のとき強制終了 
    // Neumannの安定性解析より, KAP <= 0.5でないと正しい結果が得られない
    BOOST_ASSERT(KAP <= 0.5); 

    auto x = xdiscrete();
    auto y = ydiscrete();

    auto u0 = uinitial(x, y);

    // output u(t = 0)
    if (!result_output(0, u0)){
		std::cerr << "output file not open" << std::endl;
        return -1;
    }

    for (auto t = 1; t <= TLOOP; t++){
        center_diffe(u0);

        if (!(t % (TLOOP / TREP))){
            if (!result_output(t, u0)){
		        std::cerr << "output file not open" << std::endl;
                return -1;
            }
        }
    }

    return 0;
}

namespace{
    bool result_output(int const i, my2darray const &u0){
		auto const filename = boost::format("data_nonste_diff_%d.txt") % (static_cast<int>(i / (TLOOP / TREP)));
		
        std::unique_ptr<FILE, decltype(&std::fclose)> fp(std::fopen(filename.str().c_str(), "w"), std::fclose);
        if (!fp){
            return false;
        }

        for (auto i = 0; i < XNODE; i++){
            for (auto j = 0; j < YNODE; j++){
                std::fprintf(fp.get(), "%.2f %.2f %.7f\n", static_cast<double>(i) * DX, static_cast<double>(j) * DY, u0[i][j]);
            }
            std::fprintf(fp.get(), "\n"); // for 'set pm3d' in gnuplot
        }

		return true;
    }

    myvec xdiscrete(){
        myvec x(XNODE, 0.0);

        // discrization of xrange
        for (auto i = 0; i < XNODE; i++){
            x[i] = static_cast<double>(i) * DX;
        }

        return x;
    }

    myvec ydiscrete(){
        myvec y(YNODE, 0.0);

        // discrization of yrange
        for (auto i = 0; i < YNODE; i++){
            y[i] = static_cast<double>(i) * DY;
        }

        return y;
    }

    void center_diffe(my2darray &u0){
        my2darray u(boost::extents[XNODE][YNODE]);

        switch (BCT){
            case boundary_conditon_type::DIRICLET:
                // y = 0
                for (auto i = 0; i <= XMESH; i++){
                    u[i][0] = 0.0;
                }
                // y = YRANGE
                for (auto i = 0; i <= XMESH; i++){
                    u[i][YMESH] = 0.0;
                }
                // x = 0
                for (auto i = 0; i <= YMESH; i++){
                    u[0][i] = 0.0;
                }
                // x = XRANGE
                for (auto i = 0; i <= YMESH; i++){
                    u[XMESH][i] = 0.0;
                }
                break;

            default:
                BOOST_ASSERT(!"boundary_conditon_type is abnormal value");
                break;
        }
        
        // center differential method
        for (auto i = 1; i < XNODE - 1; i++){
            for (auto j = 1; j < YNODE - 1; j++){
                u[i][j] = u0[i][j] + XKAP * (u0[i + 1][j] - 2 * u0[i][j] + u0[i - 1][j]) + YKAP * (u0[i][j + 1] - 2 * u0[i][j] + u0[i][j - 1]);
            }
        }
        
        // update value of u
        u0 = std::move(u);
    }

    my2darray uinitial(myvec const &x, myvec const &y){
        using namespace boost::math::constants;
        
        my2darray u0(boost::extents[XNODE][YNODE]);

        // initial condition of u
        for (auto i = 0; i < XNODE; i++){
            for (auto j = 0; j < YNODE; j++){
                u0[i][j] = std::sin(pi<double>() * x[i]) * std::sin(pi<double>() * y[j]);
            }
        }

        return u0;
    }
}
