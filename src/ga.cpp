#ifndef GA_CPP
#define GA_CPP

#include "boost\lexical_cast.hpp"

#include "ga.h"
#include "vector.h"
#include "support.h"

using namespace fasim;
// single-objective ===============================================================================

template<int N, int N_obj> template<typename F> void GA<N, N_obj>::run
    (F &f, Individual _lower_boundary, Individual _upper_boundary, Individual *initial_population, int initial_population_size)
{
    if(N_obj != 1)
    {
        std::cout << "GA::run: attempt to call single-objective solver with " << N_obj << " objectives\n";
        wait_and_exit();
    }

    if(!first_run)
    {
        memory_clear();
        memory_allocate();
    }
    first_run = false;

    lower_boundary = _lower_boundary;
    upper_boundary = _upper_boundary;

    seed_population(initial_population, initial_population_size);

    int
        n_crossover_children = floor(0.5 + options.crossover_fraction * (options.population_size - options.n_elite)),
        n_mutation_children = options.population_size - options.n_elite - n_crossover_children,
        n_parents = n_mutation_children + 2 * n_crossover_children;
    
    generation = 0;
    last_improvement_generation = 0;
    
    index_comparator<Objective> comp;

    for(int gen = 0; gen < options.max_generations; gen++)
    {
        // score the population
        for(int i = 0; i < options.population_size; i++)
            score[i] = f(population[i]);

        // sort score index array according to score values in ascending order
        for(int i = 0; i < options.population_size; i++)
            score_index[i] = i;

        comp.set_objective(score);
        std::sort(score_index, score_index + options.population_size, comp);

        best_score[generation] = score[score_index[0]];
        best_individual[generation] = population[score_index[0]];

        if(options.verbose)
        {
            std::cout << generation << "\t" << best_score[generation] << "\t" << population[score_index[0]] << "\n";
        }

        if(generation == options.max_generations - 1)
            break;

        if(generation > 0)
        {
            if(best_score[generation] < best_score[generation - 1])
                last_improvement_generation = generation;
        }

        // break if no improvement for last 'options.stall_generations_limit' generations


        if(generation - last_improvement_generation > options.stall_generations_limit || generation == options.max_generations)
            break;

        (this->*options.scaling)();
        
        // normalize fitness so that total fitness sum is 1.
        double cumulative_fitness = 0.;
        for(int i = 0; i < options.population_size; i++)
            cumulative_fitness += fitness[i];
        for(int i = 0; i < options.population_size; i++)
            fitness[i] /= cumulative_fitness;
        
        // get parents' indexes
        (this->*options.selection)(n_parents);

        // parents are indexes in the current generation
        std::random_shuffle(parents, parents + n_parents);

        // transfer elite children
        for(int i = 0; i < options.n_elite; i++)
            children[i] = population[score_index[i]];
        
        // perform genetic operators to get the rest of the children
        ma_step_changed = false;
        int i_parent = 0;
        for(int i = options.n_elite; i < options.n_elite + n_mutation_children; ++i, ++i_parent)
            (this->*options.mutation) (population[parents[i_parent]], children[i]);

        for(int i = options.n_elite + n_mutation_children; i < options.population_size; ++i, i_parent+=2)
            (this->*options.crossover) (population[parents[i_parent]], population[parents[i_parent + 1]], children[i]);

        for(int i = 0; i < options.population_size; ++i)
            population[i] = children[i];

        generation++;
    } // main GA loop

    std::string filename_objective = options.output_directory + "single_objective" + ".txt";
    output(filename_objective);
}


// multi-objective ================================================================================

template<int N, int N_obj> template<typename F> void GA<N, N_obj>::run_multiobjective
    (F &f, Individual _lower_boundary, Individual _upper_boundary, Individual *initial_population, int initial_population_size)
{
    if(N_obj < 2)
    {
        std::cout << "GA::run: attempt to call multi-objective solver with " << N_obj << " objective(s)\n";
        wait_and_exit();
    }

    if(first_run)
    {
        first_run = false;
    }
    else
    {
        memory_clear();
        memory_allocate();
    }

    lower_boundary = _lower_boundary;
    upper_boundary = _upper_boundary;


    int
        n_crossover_children = round(options.crossover_fraction * options.population_size),
        n_mutation_children = options.population_size - n_crossover_children,
        n_parents = n_mutation_children + 2 * n_crossover_children;

    // form the initial population and score it
    int mobj_pop_size = 2 * options.population_size;    
    seed_population(initial_population, initial_population_size);
    for(int i = 0; i < mobj_pop_size; i++)
        score[i] = f(population[i]);

    generation = 0;

    index_comparator<> comp;

    for(int gen = 0; gen < options.max_generations; gen++)
    {       
        // separate Pareto fronts =================================================================
        // unranked individuals have zero rank
        // mark all individual as unranked
        for(int i = 0; i < mobj_pop_size; i++)
            rank[i] = 0;

        // set current front rank
        total_front_ranks = 0;
        bool found_unranked = true;

        while(found_unranked)
        {
            total_front_ranks++;  
            found_unranked = false;

            // find unranked indivuduals
            for(int i = 0; i < mobj_pop_size; i++)
            {
                if(rank[i] == 0)
                {
                    found_unranked = true;
                    // test rank them with the current front rank
                    rank[i] = total_front_ranks;

                    // check all other individuals of the current front
                    for(int j = 0; j < mobj_pop_size; j++)
                    {
                        if(rank[j] == total_front_ranks && j != i)
                        {
                            // if candidate individual is dominated by some individual of the current front, delete the candidate
                            if(pareto_dominates(j, i))
                            {
                                rank[i] = 0;
                                break;
                            }
                            else
                            {
                                // if candidate individual dominates some individual of the current front, delete the latter
                                if(pareto_dominates(i, j))
                                    rank[j] = 0;
                            }
                        }
                    }

                }
            }
            
        } // while unranked individuals are found
        total_front_ranks--;

        // calculate the distances for each front =================================================
        double infinity = std::numeric_limits<double>::max();   // infinite distance flag

        for(int r = 1; r <= total_front_ranks; r++)
        {
            // get the size and indexes of the current front 
            front_size[r] = 0;
            for(int i = 0; i < mobj_pop_size; i++)
            {
                if(rank[i] == r)
                {
                    // set the distances of individuals on the front to zero
                    distance[i] = 0.;
                    // store the indexes of the individuals of the current front
                    index[front_size[r]] = i;
                    front_size[r]++;
                }
            }
            
            // accumulate distance for all objectives
            for(int obj = 0; obj < N_obj; obj++)
            {
                // get score array for current objective
                for(int i = 0; i < front_size[r]; i++)
                {
                    score_for_objective[index[i]] = score[index[i]][obj];
                }

                // sort indexes according to increasing score_for_objective_order
                comp.set_objective(score_for_objective);
                std::sort(index, index + front_size[r], comp);
                
                // distance is infinite for individuals at Pareto front ends
                distance[index[0]] = infinity;
                distance[index[front_size[r] - 1]] = infinity;

                double score_range_for_objective = fabs(score[index[0]][obj] - score[index[front_size[r] - 1]][obj]); // this can be zero, test it ==================================
                if(score_range_for_objective == 0. || score_range_for_objective == -1. * 0.)
                    score_range_for_objective = 1. + std::max(fabs(score[index[0]][obj]), fabs(score[index[front_size[r] - 1]][obj]));
                
                for(int i = 1; i < front_size[r] - 1; i++)
                {
                    if(distance[index[i]] != infinity)
                        distance[index[i]] += fabs(score[index[i + 1]][obj] - score[index[i - 1]][obj]) / score_range_for_objective;
                }
            } // for objectives ...

        } // for ranks ...

        // fill the archive =======================================================================
        int archive_size = 0;
        for(int r = 1; r <= total_front_ranks; r++)
        {
            if(front_size[r] > options.population_size - archive_size)
            {
                // archive has no room for all the front individuals, pick only the most distant
                int n_selected = options.population_size - archive_size;
                
                // get the indexes of the individuals in the current front
                int current_front_size = 0;
                for(int i = 0; i < mobj_pop_size; i++)
                {
                    if(rank[i] == r)
                    {
                        index[current_front_size] = i;
                        current_front_size++;
                    }
                }
                
                // set distance array as score
                for(int i = 0; i < front_size[r]; i++)
                {
                    score_for_objective[index[i]] = distance[index[i]];
                }

                // sort indexes according to increasing distance
                comp.set_objective(score_for_objective);
                std::sort(index, index + front_size[r], comp);

                // pick the most distant 'n_selected' individuals
                for(int i = 0; i < n_selected; i++)
                {
                    archive[archive_size] = index[front_size[r] - 1 - i];
                    archive_size++;
                }
                
                // parents selection done, exit the loop
                break;
            }
            else
            {
                // we can transfer the whole front to parents
                for(int i = 0; i < mobj_pop_size; i++)
                {
                    if(rank[i] == r)
                    {
                        archive[archive_size] = i;
                        archive_size++;
                    }
                }
            }
        } // for ranks ...

        for(int i = 0; i < options.population_size; i++)
        {
            archive_individuals[i] = population[archive[i]];
            archive_score[i] = score[archive[i]];
        }

        // calculate monitoring data ==============================================================
        // calculate average distances for non-dominated Pareto front
        average_distance[generation] = 0.;
        average_distance_deviation[generation] = 0.;

        int n_finite = 0;
        for(int i = 0; i < mobj_pop_size; i++)
        {
            if(rank[i] == 1 && distance[i] != infinity)
            {
                average_distance[generation] += distance[i];
                index[n_finite] = i;
                n_finite++;
            }
        }
        
        if(n_finite != 0)
        {
            average_distance[generation] /= n_finite;
            for(int i = 0; i < n_finite; i++)
            {
                average_distance_deviation[generation] += pow(average_distance[generation] - distance[index[i]], 2);
            }
            average_distance_deviation[generation] = sqrt(average_distance_deviation[generation] / n_finite);
        }

        // extreme Pareto solutions
        for(int obj = 0; obj < N_obj; obj++)
        {
            int i_min_elem = 0;
            for(int i = 1; i < mobj_pop_size; i++)
            {
                if(score[i][obj] < score[i_min_elem][obj])
                    i_min_elem = i;
            }

            extreme_Pareto_solution[generation][obj] = score[i_min_elem];
        }
         
        double extreme_Pareto_distance = 0.;
        if(generation > 0)
        {
            for(int i = 0; i < N_obj; i++)
                extreme_Pareto_distance += scalar_norm(extreme_Pareto_solution[generation][i] - extreme_Pareto_solution[generation - 1][i]);
        }

        if(extreme_Pareto_distance > 0 || average_distance_deviation[generation] > 0)
            spread[generation] = (extreme_Pareto_distance + average_distance_deviation[generation]) / (extreme_Pareto_distance + N_obj * average_distance[generation]);
        else
            spread[generation] = extreme_Pareto_distance;

        // output monitoring data
        if(options.verbose)
        {
            std::cout << generation << "\tspread = " << spread[generation] << "\taverage distance = " << average_distance[generation]
                << "\tPareto front size = " << front_size[1] << "\n";
        }
        
        // check termination criteria =============================================================
        bool terminate = false;
        int window = options.stall_generations_limit + 1;

        if(generation > window)
        {
            double mean_spread = 0.;
            for(int i = 0; i < window; i++)
                mean_spread += spread[generation - i];
            mean_spread /= window;

            double spread_weight = 0.5;
            double spread_change = 0.;    
            for(int i = 1; i < window; i++)
            {
                spread_change += pow(spread_weight, i) * fabs(spread[generation - i + 1] - spread[generation - i]) / (1 + spread[generation - i]);
            }
            spread_change /= window;

            if(spread_change < options.spread_change_tolerance && mean_spread >= spread[generation])
            {
                std::cout << "spread change = " << spread_change << " is less than termination tolerance\n";
                terminate = true;
            }
        }
        
        // output the data ========================================================================
        if(generation == options.max_generations - 1 || generation == 0 || generation % options.output_generations_step == 0 || terminate == true)
        {
            std::string filename_objective = options.output_directory + "objective" + boost::lexical_cast<std::string>(generation) + ".txt";
            output(filename_objective, true);

            std::string filename_parameter = options.output_directory + "parameter" + boost::lexical_cast<std::string>(generation) + ".txt";
            output(filename_parameter, false);
        }
        if(generation == options.max_generations - 1 || terminate == true)
        {
            std::string filename_objective = options.output_directory + "objective_final" + ".txt";
            output(filename_objective, true);

            std::string filename_parameter = options.output_directory + "parameter_final" + ".txt";
            output(filename_parameter, false);
        }

        
        // no need to breed in the last generation
        if(generation == options.max_generations - 1 || terminate == true)
            break;

        // breed the new generation using archive =================================================
        // get parents' indexes
        selection_tournament_multiobjective(n_parents);
        std::random_shuffle(parents, parents + n_parents);
        
        // perform genetic operators to get the children
        ma_step_changed = false; // lower the flag for adaptive mutation
        
        int i_parent = 0;
        for(int i = 0; i < n_mutation_children; ++i, ++i_parent)
            (this->*options.mutation) (population[parents[i_parent]], children[i]);

        //if(options.verbose)
        //{
        //    std::cout << "\tma_step_size = " << ma_step_size << "\n";
        //}

        for(int i = n_mutation_children; i < options.population_size; ++i, i_parent+=2)
            (this->*options.crossover) (population[parents[i_parent]], population[parents[i_parent + 1]], children[i]);

        // merge children and archive into the new population
        for(int i = 0; i < options.population_size; i++)
        {
            population[i] = children[i];
            score[i] = f(population[i]);
            population[i + options.population_size] = archive_individuals[i];
            score[i + options.population_size] = archive_score[i];
        }
        
        generation++;
    } // main GA loop
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<int N, int N_obj> bool GA<N, N_obj>::pareto_dominates(int i, int j)
{
    bool res = false;

    for(int obj = 0; obj < N_obj; obj++)
    {
        if(score[i][obj] < score[j][obj])
            res = true;         // i may dominate j (we need better value in at least one objective to dominate)

        if(score[j][obj] < score[i][obj])
            return false;       // i does not dominate j for sure
    }

    return res;
}


// scaling ========================================================================================

template <int N, int N_obj> void GA<N, N_obj>::scaling_rank()
{
    for(int i = 0; i < options.population_size; i++)
        fitness[score_index[i]] = 1. / sqrt(i + 1.);  
}


// selection ======================================================================================

template <int N, int N_obj> void GA<N, N_obj>::selection_stochastic_uniform(int n)
{
    double *wheel = new double[options.population_size];
    
    wheel[0] = fitness[0];
    for(int i = 1; i < options.population_size; ++i)
        wheel[i] = wheel[i - 1] +  fitness[i];

    double step = 1. / (double)n,
           position = dist01(rnd_generator) * step;
    int j_low = 0;
    
    for(int i = 0; i < n; ++i)
    {
        for(int j = j_low; j < options.population_size; ++j)
        {
            if(position < wheel[j])
            {
                parents[i] = j;
                j_low = j;
                break;
            }
        }
        position += step;
    }
}


template <int N, int N_obj> void GA<N, N_obj>::selection_tournament_shuffle(int n)
{
    int *indexes = new int [options.population_size];
    for(int i = 0; i < options.population_size; ++i)
        indexes[i] = i;

    for(int i = 0; i < n; i++)
    {
        std::random_shuffle(indexes, indexes + options.population_size);
        
        parents[i] = indexes[0];
        for(int j = 1; j < options.tournament_size; j++)
        {
            if(fitness[indexes[j]] > fitness[parents[i]])
                parents[i] = indexes[j];
        }
    }

    delete [] indexes;
}


template <int N, int N_obj> void GA<N, N_obj>::selection_tournament(int n)
{
    for(int i = 0; i < n; i++)
    {
        int best = floor(options.population_size * dist01(rnd_generator));

        for(int j = 1; j < options.tournament_size; j++)
        {
            int candidate = floor(options.population_size * dist01(rnd_generator));
            
            if(fitness[candidate] > fitness[best])
                best = candidate;
        }
        parents[i] = best;
    }
}


template <int N, int N_obj> void GA<N, N_obj>::selection_tournament_multiobjective(int n)
{
    for(int i = 0; i < n; i++)
    {
        int best = floor(options.population_size * dist01(rnd_generator));
        best = archive[best];

        for(int j = 1; j < options.tournament_size_multiobjective; j++)
        {
            int candidate = floor(options.population_size * dist01(rnd_generator));
            candidate = archive[candidate];
            
            if(rank[candidate] < rank[best])
                best = candidate;
            else
            {
                if(rank[candidate] == rank[best] && distance[candidate] > distance[best])
                    best = candidate;
            }
        }
        parents[i] = best;
    }
}

// crossover ======================================================================================

template <int N, int N_obj> void GA<N, N_obj>::crossover_arithmetic
	(const Individual &parent1, const Individual &parent2, Individual &child)
{
    double R = dist01(rnd_generator);
    child = R * parent1 + (1. - R) * parent2;
}


template <int N, int N_obj> void GA<N, N_obj>::crossover_scattered
	(const Individual &parent1, const Individual &parent2, Individual &child)
{
    for(int i = 0; i < N; ++i)
    {
        if(dist01(rnd_generator) >= 0.5)
        {
            child[i] = parent1[i];
        }
        else
        {
            child[i] = parent2[i];
        }
    }
    
    if( !feasible(child) )
        crossover_arithmetic(parent1, parent2, child);
}


template <int N, int N_obj> void GA<N, N_obj>::crossover_BLX
	(const Individual &parent1, const Individual &parent2, Individual &child)
{
    for(int i = 0; i < N; i++)
    {
        double R = dist01(rnd_generator);

        double I = parent2[i] - parent1[i],
               low = parent1[i] - I * options.crossover_BLX_alpha,
               high = parent2[i] + I * options.crossover_BLX_alpha;

        child[i] = R * low + (1. - R) * high;

        // modify child to satisfy constraints
        //child[i] = std::max(child[i], lower_boundary[i]);
        //child[i] = std::min(child[i], upper_boundary[i]);
        if(child[i] < lower_boundary[i] || child[i] > upper_boundary[i])
            child[i] = (R * parent1[i] + (1. - R) * parent2[i]);
    }
}


// mutation =======================================================================================

template <int N, int N_obj> void GA<N, N_obj>::mutation_gaussian
	(const Individual &parent, Individual &child)
{
    // this mutation breaks constraints
    double scale_factor = options.mutation_gaussian_scale * (1. - options.mutation_gaussian_shrink * generation / options.max_generations);
    Individual scale = scale_factor * (upper_boundary - lower_boundary);

    for(int i = 0; i < N; ++i)
        child[i] = parent[i] + scale[i] * normal01(rnd_generator);
}


template <int N, int N_obj> void GA<N, N_obj>::mutation_adaptive
	(const Individual &parent, Individual &child)
{
    double tol = 1.e-8;             // tolerance, pure magic

    // set step size
    if(generation <= 1)
    {
        ma_step_size = 1.;
    }
    else
    {
        if(!ma_step_changed)
        {
            ma_step_changed = true;

            if(N_obj == 1)
            {
                if(best_score[generation] < best_score[generation - 1])
                    ma_step_size = std::min(1., ma_step_size * 4.);
                else
                    ma_step_size = std::max(tol, ma_step_size / 4.);
            }
            else
            {
                if(spread[generation] > spread[generation - 1])
                    ma_step_size = std::min(1., ma_step_size * 4.);
                else
                    ma_step_size = std::max(tol, ma_step_size / 4.);
            }
            
        }
    }

    // set logarithmic scale
    Individual scale;
    for(int i = 0; i < N; i++)
    {
        double exponent = 0.;
        if( fabs(lower_boundary[i]) > tol )
            exponent += 0.5 * log(fabs(lower_boundary[i])) / log(2.);
        if( fabs(upper_boundary[i]) > tol )
            exponent += 0.5 * log(fabs(upper_boundary[i])) / log(2.);

        scale[i] = pow(2., exponent);
    }

    Individual raw_basis[N],
               basis[N],
               tangent_cone[N],
               dir[2 * N];
    double dir_sign[4 * N];

    int n_tangent = 0,
        n_basis = N,
        index_vector[4 * N],
        order_vector[4 * N];

    // calculate mutation direction set
    // tangent components
    for(int i = 0; i < N; i++)
    {
        if( fabs(parent[i] - lower_boundary[i]) < tol || fabs(parent[i] - upper_boundary[i]) < tol )
        {
            tangent_cone[n_tangent] = 0.;
            tangent_cone[n_tangent][i] = 1.;
            n_tangent++;
        }
    }
    
    // raw basis vectors
    double poll_param = 1. / sqrt(ma_step_size);

    for(int i = 0; i < n_basis; i++)
    {
        raw_basis[i] = 0.;
        raw_basis[i][i] = poll_param * (dist01(rnd_generator) >= 0.5 ? 1. : -1.);
        for(int j = i + 1; j < n_basis; j++)
        {
            raw_basis[i][j] = round( (poll_param + 1.) * dist01(rnd_generator) - 0.5);
        }
    }

    // basis as random permutation of raw basis
    for(int j = 0; j < n_basis; j++)
        order_vector[j] = j;
    std::random_shuffle(order_vector, order_vector + n_basis);

    for(int i = 0; i < n_basis; i++)
        for(int j = 0; j < n_basis; j++)
            basis[i][j] = raw_basis[order_vector[i]][order_vector[j]];

    // prerare random direction mutation
    int n_dir = n_tangent + n_basis;

    for(int i = 0; i < n_basis; i++)
        dir[i] = basis[i];

    for(int i = 0; i < n_tangent; i++)
        dir[n_basis + i] = tangent_cone[i];

    for(int i = 0; i < n_basis; i++)
    {
        index_vector[i] = i;
        dir_sign[i] = 1.;
    }
    
    int i_base = n_basis;
    for(int i = i_base; i < i_base + n_basis; i++)
    {
        index_vector[i] = i - i_base;
        dir_sign[i] = -1.;
    }
    
    i_base += n_basis;
    for(int i = i_base; i < i_base + n_tangent; i++)
    {
        index_vector[i] = i - i_base;
        dir_sign[i] = 1.;
    }
    
    i_base += n_tangent;
    for(int i = i_base; i < i_base + n_tangent; i++)
    {
        index_vector[i] = i - i_base;
        dir_sign[i] = -1.;
    }
    
    int n_dir_total = 2 * n_dir;
    for(int i = 0; i < n_dir_total; i++)
        order_vector[i] = i;
    std::random_shuffle(order_vector, order_vector + n_dir_total);

    // finally, mutate
    double success = false;
    for(int i = 0; i < n_dir_total; i++)
    {
        int k = index_vector[order_vector[i]];
        Individual direction = dir_sign[k] * dir[k];
        
        child = parent + ma_step_size * scale * direction;

        if(feasible(child))
        {
            success = true;
            break;
        }
    }

    if( !success )
    {
        child = parent;
        if(options.verbose && N < 5)
            std::cout << "mutation failed at x = " << parent << "\n";
    }
}


// all the rest ===================================================================================

template<int N, int N_obj> typename GA<N, N_obj>::Individual GA<N, N_obj>::random_individual
    (const typename Individual &lower_boundary, const typename Individual &upper_boundary)
{
    Individual res;
    for(int i = 0; i < N; ++i)
        res[i] = lower_boundary[i] + dist01(rnd_generator) * (upper_boundary[i] - lower_boundary[i]);
    
    return res;
}


template<int N, int N_obj>  void GA<N, N_obj>::seed_population(Individual *initial_population = NULL, int initial_population_size = 0)
{
    // seed double the size for multiobjective optimization
    double size = (N_obj > 1 ? 2 * options.population_size : options.population_size);
    
    if(initial_population == NULL)
    {
        for(int i = 0; i < size; ++i)
            population[i] = random_individual(lower_boundary, upper_boundary);
    }
    else
    {
        int initial_size = (initial_population_size == 0 ? options.population_size: initial_population_size);
 
        for(int i = 0; i < initial_size; i++)
        {
            population[i] = initial_population[i];
            if( !feasible(population[i]) )
            {
                std::cout << "GA::run: initial population individual " << i << " is violating boundaries\n";
                wait_and_exit();
            }
        }

        for(int i = initial_size; i < size; i++)
            population[i] = random_individual(lower_boundary, upper_boundary);
    }
}


template<int N, int N_obj> bool GA<N, N_obj>::feasible(const Individual &x)
{
    for(int i = 0; i < N; i++)
    {
        if(x[i] > upper_boundary[i] || x[i] < lower_boundary[i])
            return false;
    }

    return true;
}


template <int N, int N_obj> void GA<N, N_obj>::output(const std::string &filename, bool output_objective)
{
    std::ofstream file_output(filename.c_str());
		
	if (!file_output) 
	{
		std::cout << "f_write_table: error reading " << filename.c_str() << "\n";
		wait_and_exit();
	}

    if(N_obj == 1)
    {
        for(int i = 0; i <= generation; i++)
        {
            file_output << best_score[i] << "\t" << best_individual[i];
            if(i != generation)
                file_output << "\n";
        }
    }
    else
    {
        bool first = true;
        int n_out = 2. * options.population_size;	
        for(int i = 0; i < n_out; i++)
	    {
            if(rank[i] == 1)
            {
                if(first)
                {
                    if(output_objective)
                        file_output << score[i];
                    else
                        file_output << population[i];
                    first = false;
                }
                else
                {
                    file_output << "\n";
                    if(output_objective)
                        file_output << score[i];
                    else
                        file_output << population[i];                  
                }
            }
	    }
    }
	file_output.close();
}


#endif