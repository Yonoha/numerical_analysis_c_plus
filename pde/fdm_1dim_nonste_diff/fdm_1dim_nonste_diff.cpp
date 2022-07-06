// 1-dim nonsteady diffusion equation by fdm

#include <cmath>    // for std::sin
#include <fstream>  // for std::ofstream 
#include <iostream> // for std::cerr, std::endl
#include <vector>   // for std::vector
#include <boost/assert.hpp>                   // for BOOST_ASSERT
#include <boost/math/constants/constants.hpp> // for boost::math::constants::pi

namespace{
    using myvec = std::vector<double>;

    static auto constexpr TEND = 10.0;
    static auto constexpr TLOOP = 10000;
    static auto constexpr DT = TEND / TLOOP;
    static auto constexpr TREP = TLOOP / 1000;

    static auto constexpr MESH = 100;
    static auto constexpr XRANGE = 1.0;
    static auto constexpr DX = XRANGE / MESH;

    // boundary condition
    static auto constexpr BOUND_LEFT = 0.0;
    static auto constexpr BOUND_RIGHT = 0.0;
    
    // coefficient
    static auto constexpr DIFF = 0.01;
    static auto constexpr KAP = DIFF * DT / (DX * DX);
    
    // result output file
    bool output_result(std::vector<myvec> const &uxt);

    // discretization of x
    myvec discrete();

    // initial condition u(x, 0)
    myvec initial(myvec const &x);

    // center differential method
    std::vector<myvec> differencial_method(myvec const &u0);
}

int main(){
    // KAP > 0.5のとき強制終了 
    // Neumannの安定性解析より, KAP <= 0.5でないと正しい結果が得られない
    BOOST_ASSERT(KAP <= 0.5);

    auto x = discrete();

    auto u0 = initial(x);

    auto uxt = differencial_method(u0);

    output_result(uxt);

    if (!output_result(uxt)){
        std::cerr << "error! output file not open" << std::endl;
        return -1;
    }

    return 0;
}

namespace{
    bool output_result(std::vector<myvec> const &uxt){
        std::vector<myvec> u(uxt.size(), (myvec(TREP + 1)));
        
        // output only TREP + 1 (TLOOP = 0, 1000, ... , 10000)
        for (auto i = 0; i <= TREP; i++){
            for (auto j = 0; j <= MESH; j++){
                u[j][i] = uxt[j][1000 * i];
            }
        }

        std::ofstream ofs;
        ofs.open("data_center_diffe_diff.txt");

        if (!ofs){
            return false;
        }

        for (auto j = 0; j <= MESH; j++){
            ofs << static_cast<double>(j) * DX << " ";
            
            for (auto i = 0; i <= TREP; i++){
                ofs << u[j][i] << " ";
            }
            ofs << std::endl;
        }

        return true;
    }

    myvec discrete(){
        myvec x(MESH + 1, 0.0);

        // discrization of xrange
        for (int i = 0; i <= MESH; i++){
            x[i] = i * DX;
        }

        return x;
    }

    myvec initial(myvec const &x){
        myvec u0(MESH + 1, 0.0);

        // initial condition
        for (int i = 0; i<= MESH; i++){
            u0[i] = std::sin(boost::math::constants::pi<double>() * x[i]);
        }

        return u0;
    }

    std::vector<myvec> differencial_method(myvec const &u0){
        std::vector<myvec> uxt(MESH + 1, (myvec(TLOOP + 1)));

        // initial condition
        for (auto j = 0; j <= MESH; j++){
            uxt[j][0] = u0[j];
        }

        for (auto i = 0; i < TLOOP; i++){
            // boundary condition
            uxt[0][i + 1] = BOUND_LEFT;
            uxt[MESH][i + 1] = BOUND_RIGHT;

            // center differential method
            for (auto j = 1; j < MESH; j++){
                uxt[j][i + 1] = uxt[j][i] + KAP * (uxt[j + 1][i] - 2.0 * uxt[j][i] + uxt[j - 1][i]);
            }
        }

        return uxt;
    }
}
