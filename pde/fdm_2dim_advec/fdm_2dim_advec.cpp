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
    using mypair = std::tuple<int, int>;
    using my2darray = boost::multi_array<double, 2>;
    using myvec = std::vector<double>;

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
    static auto constexpr TLOOP = 10000; // (TLOOP / TREP) loopごとにgifにするのでTLOOP > 10 * TREP程度が良い
    static auto constexpr DT = TEND / TLOOP;

    static auto constexpr XRANGE = 1.0;
    static auto constexpr XMESH = 100;
    static auto constexpr XNODE = XMESH + 1;
    static auto constexpr DX = XRANGE / XMESH; 

    static auto constexpr YRANGE = 1.0;
    static auto constexpr YMESH = 100;
    static auto constexpr YNODE = YMESH + 1;
    static auto constexpr DY = YRANGE / YMESH;

    static auto constexpr XADVEC = 1.0; // x direction advection number
    static auto constexpr YADVEC = -1.0; // y direction advection number
    static auto constexpr XCOURANT = XADVEC * DT / DX;
    static auto constexpr YCOURANT = YADVEC * DT / DY;
    static auto constexpr COURANT = XCOURANT + YCOURANT;

    // output result for each i
    bool result_output(int const t, my2darray const &u0);

    // initial condition of u(x, y, t)
    my2darray uinitial(myvec const &x, myvec const &y);

    // make down and up value of x[i] for upwind
    mypair make_param_du(int const j);

    // make left and right value of x[i] for upwind
    mypair make_param_lr(int const i);

    // discrezation of x and y
    myvec xdiscrete();
    myvec ydiscrete();

    // 1st accuracy upwind differnece method
    void upwind_diffe(my2darray &u0);
}

int main(){
    // std::fabs(COURANT) > 1.0のとき強制終了 
    // CFL条件より, std::fabs(COURANT) <= 1.0でないと正しい結果が得られない
    BOOST_ASSERT(std::fabs(COURANT) <= 1.0);

    auto x = xdiscrete();
    auto y = ydiscrete();

    auto u0 = uinitial(x, y);

    // output u(t = 0)
    if (!result_output(0, u0)){
		std::cerr << "output file not open" << std::endl;
        return -1;
    }

    for (auto t = 1; t <= TLOOP; t++){
        upwind_diffe(u0);

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
    bool result_output(int const t, my2darray const &u0){
		auto const filename = boost::format("data_fdm_2dim_nonste_advec_%d.txt") % (static_cast<int>(t / (TLOOP / TREP)));
		
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

    mypair make_param_du(int const j){
        auto down = j - 1;
        auto up = j + 1;

        // periodic condition
        if (j == 0){
            down = YMESH;
        }

        else if (j == YMESH){
            up = 0;
        } 

        return std::make_pair(down, up);
    }

    mypair make_param_lr(int const i){
        auto left = i - 1;
        auto right = i + 1;

        // periodic condition
        if (i == 0){
            left = XMESH;
        }

        else if (i == XMESH){
            right = 0;
        } 

        return std::make_pair(left, right);
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

    void upwind_diffe(my2darray &u0){
        my2darray u(boost::extents[XNODE][YNODE]);
        my2darray xdiff(boost::extents[XNODE][YNODE]);
        my2darray ydiff(boost::extents[XNODE][YNODE]);

        switch (BCT){
            case boundary_conditon_type::DIRICLET:
                break;

            case boundary_conditon_type::LEFT_NEUMANN:
                break;

            case boundary_conditon_type::RIGHT_NEUMANN:
                break;

            case boundary_conditon_type::PERIODIC:
                for (auto i = 0; i <= XMESH; i++){
                    for (auto j = 0; j <= YMESH; j++){
                        auto const [left, right] = make_param_lr(i);
                        auto const [down, up] = make_param_du(j);

                        // 1st accuracy upwind difference method
                        xdiff[i][j] = ((XCOURANT + std::fabs(XCOURANT)) / 2.0) * (u0[i][j] - u0[left][j]) + ((XCOURANT - std::fabs(XCOURANT)) / 2.0) * (-u0[i][j] + u0[right][j]);

                        ydiff[i][j] = ((YCOURANT + std::fabs(YCOURANT)) / 2.0) * (u0[i][j] - u0[i][down]) + ((YCOURANT - std::fabs(YCOURANT)) / 2.0) * (-u0[i][j] + u0[i][up]);

                        u[i][j] = u0[i][j] - xdiff[i][j] - ydiff[i][j];
                    }
                }
                break;

            default:
                BOOST_ASSERT(!"boundary_conditon_type is abnormal value");
                break;
        }
        
        // update value of u
        u0 = std::move(u);
    }
}
