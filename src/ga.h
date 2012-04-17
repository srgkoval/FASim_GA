#ifndef GA_H
#define GA_H

#include <iostream>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost\random\variate_generator.hpp>

#include <boost/math/constants/constants.hpp>

#include "vector.h"

using namespace fasim;

// forward declarations
template <int N, int N_obj> class GA;
template <int N, int N_obj> class GA_options;  


template <int N, int N_obj> class GA_options
{
public:
    bool multiobjective_problem;

    int population_size,
        n_elite;

    int max_generations,
        stall_generations_limit;

    int tournament_size,
        tournament_size_multiobjective;

    double
        crossover_fraction,
        crossover_BLX_alpha,
        mutation_gaussian_scale,
        mutation_gaussian_shrink;

    double pareto_fraction,
           pareto_R;

    double spread_change_tolerance;     // termination criterion

    typename GA<N, N_obj>::pFitnessScaling scaling;
    typename GA<N, N_obj>::pSelection selection;
    typename GA<N, N_obj>::pCrossover crossover;
    typename GA<N, N_obj>::pMutation mutation;

    bool verbose;
    std::string output_directory;
    int output_generations_step;

    GA_options()
    {
        if(N_obj == 1)
        {
            population_size = 30;
            max_generations = 100;
        }
        else
        {
            population_size = 50;
            max_generations = 250;
        }

        n_elite = 2;

        stall_generations_limit = 50;

        tournament_size = 2;
        tournament_size_multiobjective = 2;

        crossover_fraction = 0.8;

        crossover_BLX_alpha = 0.5;
        mutation_gaussian_scale = 0.5;
        mutation_gaussian_shrink = 0.75;

        pareto_fraction = 0.35;
        pareto_R = 0.8;

        spread_change_tolerance = 1.e-4;

        scaling = &GA<N, N_obj>::scaling_rank;
        selection = &GA<N, N_obj>::selection_tournament;
        crossover = &GA<N, N_obj>::crossover_BLX;
        mutation = &GA<N, N_obj>::mutation_adaptive;

        verbose = true;
        output_directory = "e:\\Programming - code and etc\\Genetic Algorithm\\Output\\";
        output_generations_step = 5;
    }
};


// comparator class to emulate MATLAB's [~,i] = sort(scores) functionality
// allows sorting of int index array according to score values that correspond to these indexes
// no memory management needed in this class, all allocations/deallocations are handled elsewhere
template <typename T = double> class index_comparator
{
public:
    T *score;
    index_comparator(): score(NULL) {}

    void set_objective(T *_score) {score = _score;}
    
    bool operator()(int l, int r)
    {
        return score[l] < score[r];
    }
};


template <int N, int N_obj = 1> class GA
{
public:
    typedef Vector<double, N> Individual;
    typedef Vector<double, N_obj> Objective;

    typedef void FitnessScaling();
    typedef void(GA<N, N_obj>::*pFitnessScaling)();

    typedef void Selection(int n);
    typedef void (GA<N, N_obj>::*pSelection)(int n);

    typedef void Crossover(const Individual &parent1, const Individual &parent2, Individual &child);
    typedef void (GA<N, N_obj>::*pCrossover)(const Individual &parent1, const Individual &parent2, Individual &child);

    typedef void Mutation(const Individual &parent, Individual &child);
    typedef void (GA::*pMutation)(const Individual &parent, Individual &child);
    
    GA_options<N, N_obj> options;

    int generation,
        last_improvement_generation;

    Individual *population,
               *children,
               *best_individual;
    int *parents;
    
    Objective *score,
              *best_score;
    
    double *fitness;

    int *score_index;

    // multiobjective data
    int *archive;
    Individual *archive_individuals;
    Objective *archive_score;

    int total_front_ranks;
    int *front_size;
    int *rank;
    double *distance;

    bool pareto_dominates(int i, int j);   // check if i-th individual Pareto dominates j-th    

    int *index;                    // temporary
    double *score_for_objective;        // temporary

    double *average_distance,
           *average_distance_deviation,
           *spread;

    Vector<Objective, N_obj> *extreme_Pareto_solution;

    // ---


    Individual lower_boundary,
               upper_boundary;

    FitnessScaling scaling_rank;

    Selection selection_stochastic_uniform,
              selection_tournament,
              selection_tournament_shuffle,
              selection_tournament_multiobjective;

    Crossover crossover_arithmetic,
              crossover_scattered,
              crossover_BLX;

    Mutation mutation_gaussian,
             mutation_adaptive;
    
    double ma_step_size;                                    // step size for mutation adaptive
    bool ma_step_changed;                                      // flag for adaptive mutation to change step once in generation

    mutable boost::random::mt19937 rnd_generator;
	boost::random::uniform_real_distribution<> dist01;
    boost::random::normal_distribution<> normal01;
    boost::random::uniform_int_distribution<> uniform_int0N;
    boost::random::variate_generator<boost::random::mt19937, boost::random::uniform_int_distribution<>> int_gen;

    bool first_run;
    void memory_allocate();
    void memory_clear();

    GA();
    ~GA();

    Individual random_individual(const typename Individual &lower_boundary, const typename Individual &upper_boundary);
    bool feasible(const Individual &x);

    void seed_population(Individual *initial_population = NULL, int initial_population_size = 0);
    
    template<typename F> void run(F &f, Individual _lower_boundary, Individual _upper_boundary, Individual *initial_population = NULL, int initial_population_size = 0);
    template<typename F> void run_multiobjective(F &f, Individual _lower_boundary, Individual _upper_boundary,
        Individual *initial_population = NULL, int initial_population_size = 0);

    void output(const std::string &filename, bool output_objective = true);
};


template<int N, int N_obj> void GA<N, N_obj>::memory_allocate()
{
    population = new Individual [2 * options.population_size];
    children =  new Individual [options.population_size];
    best_individual =  new Individual [options.max_generations];

    score = new Objective [2 * options.population_size];
    fitness = new double [options.population_size];
    best_score = new Objective [options.max_generations];
    score_index = new int[options.population_size];
    parents = new int [2 * options.population_size];

    if(N_obj > 1)
    {
        archive = new int [options.population_size];
        archive_individuals = new Individual [options.population_size];
        archive_score = new Objective [options.population_size];
  
        front_size = new int [2 * options.population_size];
        rank = new int [2 * options.population_size];
        distance = new double [2 * options.population_size];

        index = new int [2 * options.population_size];
        score_for_objective = new double [2 * options.population_size];

        average_distance = new double [options.max_generations];
        average_distance_deviation = new double [options.max_generations];
        spread = new double [options.max_generations];
        extreme_Pareto_solution = new Vector<Objective, N_obj>  [options.max_generations];
    }
}


template<int N, int N_obj> void GA<N, N_obj>::memory_clear()
{
    delete [] population;
    delete [] children;
    delete [] best_individual;

    delete [] score;
    delete [] fitness;
    delete [] best_score;
    delete [] score_index;
    delete [] parents;

    if(N_obj > 1)
    {
        delete [] archive;
        delete [] archive_individuals;
        delete [] archive_score;

        delete [] front_size;
        delete [] rank;
        delete [] distance;

        delete [] index;
        delete [] score_for_objective;

        delete [] average_distance;
        delete [] average_distance_deviation;
        delete [] spread;
        delete [] extreme_Pareto_solution;
    }
}


template<int N, int N_obj> GA<N, N_obj>::GA()
    : dist01(0., 1.), generation(0), int_gen(rnd_generator, uniform_int0N)
{
    rnd_generator.seed(static_cast<int>(std::time(NULL)));
    srand(unsigned(time(NULL)));
        
    memory_allocate();
    first_run = true;
}


template<int N, int N_obj> GA<N, N_obj>::~GA()
{
    memory_clear();
}



#include "ga.cpp"

#endif