#ifndef _MCPELE_METROPOLIS_TEST_H__
#define _MCPELE_METROPOLIS_TEST_H__

#include <random>

#include "pele/array.h"
#include "mc.h"

namespace mcpele {

/*! \class MetropolisTest
 *  \brief Metropolis acceptance test.
 *
 *  The Metropolis acceptance criterion accepts each move with probability
 *  @f[
 *  P( x_{old} \Rightarrow x_{new}) = min \{ 1, \exp [- \beta (E_{new} - E_{old})] \}
 *  @f]
 *  where \f$\beta\f$ is the reciprocal of the temperature
 */

class MetropolisTest : public AcceptTest {
protected:
    size_t m_seed;                                          /*!< Seed for random number generator*/
    std::mt19937_64 m_generator;                            /*!< Mersenne twister random number generator*/
    std::uniform_real_distribution<double> m_distribution;  /*!< Distribution for random number generator*/
public:
    /*!Constructor.
      \param rseed a random positive integer
    */
    MetropolisTest(const size_t rseed);

    /*!Destructor.*/
    virtual ~MetropolisTest() {}

    /*!
     * \param trial_coords array of trial coordinates
     * \param trial_energy energy associated with trial_coords
     * \param old_coords array of coordinates before TakeStep
     * \param old_energy energy associeted with old_coords
     * \param temperature system temperature
     * \param mc pointer to the MC runner
     * \return bool for outcome of the test
     * */
    virtual bool test(pele::Array<double> &trial_coords, double trial_energy,
            pele::Array<double> & old_coords, double old_energy, double temperature,
            MC * mc);

    /*!\return seed used for rng*/
    size_t get_seed() const {return m_seed;}
    void set_generator_seed(const size_t inp) { m_generator.seed(inp); }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_METROPOLIS_TEST_H__
