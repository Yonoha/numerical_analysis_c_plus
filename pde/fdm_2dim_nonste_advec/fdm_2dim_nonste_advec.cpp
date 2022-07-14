// 2dim non-steady advection equation by fdm
// use upwind difference method 

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
    static auto constexpr BCT = boundary_conditon_type::PERIODIC;

    static auto constexpr TEND = 1.0;
    static auto constexpr TREP = 100;
    static auto constexpr TLOOP = 10000; // TREP loopごとにgifにするのでTLOOP > 10 * TREP程度が良い
    static auto constexpr DT = TEND / TLOOP;

    static auto constexpr XRANGE = 1.0;
    static auto constexpr XMESH = 100;
    static auto constexpr XNODE = XMESH + 1;
    static auto constexpr DX = XRANGE / XMESH; 

    static auto constexpr YRANGE = 1.0;
    static auto constexpr YMESH = 100;
    static auto constexpr YNODE = YMESH + 1;
    static auto constexpr DY = YRANGE / YMESH;

    static auto constexpr XADVEC = 0.5; // x direction advection number
    static auto constexpr YADVEC = 1.0; // y direction advection number
    static auto constexpr XCAURANT = XADVEC * DT / DX;
    static auto constexpr YCAURANT = YADVEC * DT / DY;
    static auto constexpr ADVEC = XCAURANT + YCAURANT;

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
    // ADVEC > 0.5のとき強制終了 
    // CLT条件より, ADVEC <= 1.0でないと正しい結果が得られない
    BOOST_ASSERT(0.0 <= ADVEC && ADVEC <= 1.0); 

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
		auto const filename = boost::format("data_fdm_2dim_nonste_advec_%d.txt") % (static_cast<int>(i / (TLOOP / TREP)));
		
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

    my2darray uinitial(myvec const &x, myvec const &y){
        using namespace boost::math::constants;
        
        my2darray u0(boost::extents[XNODE][YNODE]);

        // initial condition of u
        for (auto i = 0; i < XNODE; i++){
            for (auto j = 0; j < YNODE; j++){
                u0[i][j] = std::sin(2.0 * pi<double>() * x[i]) * std::sin(2.0 * pi<double>() * y[j]);
            }
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
                for (auto i = 1; i < YMESH; i++){
                    u[0][i] = 0.0;
                }

                // x = XRANGE
                for (auto i = 1; i < YMESH; i++){
                    u[XMESH][i] = 0.0;
                }
                break;

            case boundary_conditon_type::LEFT_NEUMANN:

            case boundary_conditon_type::RIGHT_NEUMANN:

            case boundary_conditon_type::PERIODIC:
                // y = 0
                u[0][0] = u0[0][0] - XCAURANT * (u0[0][0] - u0[XMESH][0]) - YCAURANT * (u0[0][0] - u0[0][YMESH]);

                for (auto i = 1; i <= XMESH; i++){
                    u[i][0] = u0[i][0] - XCAURANT * (u0[i][0] - u0[i - 1][0]) - YCAURANT * (u0[i][0] - u0[i][YMESH]);
                }

                // y = YRANGE
                u[0][YMESH] = u0[0][YMESH] - XCAURANT * (u0[0][YMESH] - u0[XMESH][YMESH]) - YCAURANT * (u0[0][YMESH] - u0[0][YMESH - 1]);

                for (auto i = 1; i <= XMESH; i++){
                    u[i][YMESH] = u0[i][YMESH] - XCAURANT * (u0[i][YMESH] - u0[i - 1][YMESH]) - YCAURANT * (u0[i][YMESH] - u0[i][YMESH - 1]);
                }

                // x = 0
                for (auto i = 1; i < YMESH; i++){
                    u[0][i] = u0[0][i] - XCAURANT * (u0[0][i] - u0[XMESH][i]) - YCAURANT * (u0[0][i] - u0[0][i - 1]);
                }

                // x = XRANGE
                for (auto i = 1; i < YMESH; i++){
                    u[XMESH][i] = u0[XMESH][i] - XCAURANT * (u0[XMESH][i] - u0[XMESH - 1][i]) - YCAURANT * (u0[XMESH][i] - u0[XMESH][i - 1]);
                }
                break;

            default:
                BOOST_ASSERT(!"boundary_conditon_type is abnormal value");
                break;
        }
        
        // center differential method
        for (auto i = 1; i < XMESH; i++){
            for (auto j = 1; j < YMESH; j++){
                u[i][j] = u0[i][j] - XCAURANT * (u0[i][j] - u0[i - 1][j]) - YCAURANT * (u0[i][j] - u0[i][j - 1]);
            }
        }
        
        // update value of u
        u0 = std::move(u);
    }
}
