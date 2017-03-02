

// ================================================================================

void Cloth::Animate() {

    //all springs in current mesh
    //a variable to trigger recaculate timestep
    int re_calculate = 0;
    //Signal 0: initial
    //Signal 1: Good to go
    //Signal 2: need recaculate
    while (re_calculate == 0 || re_calculate == 2) {
        if(re_calculate == 2){
            //reset to initial before recalculate
            re_calculate = 0;
        }
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                if (getParticle(i, j).isFixed() != true) {
                    //store all particles connected with
                    //current particle with specific type of spring
                    std::vector<std::pair<ClothParticle *, float>> structural;
                    std::vector<std::pair<ClothParticle *, float>> shear;
                    std::vector<std::pair<ClothParticle *, float>> flexion;
                    std::vector<std::tuple<ClothParticle *, float, int>> all_springs;
                    all_springs = identify_strings(i, j, structural, shear, flexion, springs);
                    glm::vec3 internal_force = glm::vec3(0, 0, 0);

                    //calculate internal force
                    for (int z = 0; z < all_springs.size(); ++z) {
                        //iterate all springs
                        glm::vec3 pkl = std::get<0>(all_springs[z])->getOriginalPosition();
                        glm::vec3 pij = getParticle(i, j).getOriginalPosition();
                        glm::vec3 tmp = pkl - pij;
                        float l0 = std::get<1>(all_springs[z]);
                        double k = 0;
                        if (std::get<2>(all_springs[z]) == 0) {
                            k = k_structural;
                        } else if (std::get<2>(all_springs[z]) == 1) {
                            k = k_shear;
                        } else if (std::get<2>(all_springs[z]) == 2) {
                            k = k_bend;
                        }
                        internal_force += (float(k) * (tmp - glm::normalize(tmp) * l0));
                    }

                    //calculate gravity
                    getParticle(i, j).setAcceleration(args->gravity);
                    glm::vec3 external_force = getParticle(i, j).getForce();

                    //calculate viscous damping
                    glm::vec3 damp_force = -(float(damping) * getParticle(i, j).getOriginalVolicity());

                    //update particles

                    glm::vec3 aceleration =
                            (internal_force + external_force + damp_force) / float(getParticle(i, j).getMass());
                    getParticle(i, j).setAcceleration(aceleration);
                    glm::vec3 velocity = getParticle(i, j).getOriginalVolicity() + float(args->timestep) * aceleration;
                    float vol = glm::length(float(args->timestep) * aceleration);

                    //If speed difference larger than certain therehold, recalculate time step
                    if(vol > 5.0){
                        re_calculate = 2;
                        args->timestep = args->timestep/2.0;
                        break;
                    }
                    getParticle(i, j).setVelocity(velocity);
                    glm::vec3 position = getParticle(i, j).getOriginalPosition() + float(args->timestep) * velocity;
                    getParticle(i, j).setPosition(position);

                    //if reach the last element and nothing exceed msc vol. we are good to go.
                    if(i == nx-1 && j == ny-1){
                        re_calculate = 1;
                    }
                }
            }
            if(re_calculate == 2){
                break;
            }
        }
    }

    for (int i=0; i<6; i++) {
        std::unordered_map<std::string, std::tuple<ClothParticle *, ClothParticle *, int, float> >::iterator it = springs.begin();
        while (it != springs.end()) {
            //adjust springs here.
            float standard_length = std::get<3>(it->second);
            float actual_length = glm::length(
                    std::get<1>(it->second)->getPosition() - std::get<0>(it->second)->getPosition());
            float provot_correction = 0;
            if (std::get<2>(it->second) == 0) {
                provot_correction = provot_structural_correction;
            } else if (std::get<2>(it->second) == 1) {
                provot_correction = provot_shear_correction;
            }
            //100 means do not do any correction
            if (provot_correction == 100) {
                it++;
                continue;
            }

            float max_length = standard_length * (1.0f + provot_correction);
            //test if need to be corrected
            if (actual_length < max_length) {
                it++;
                continue;
            }

            //working on here
            float correction_factor_x = max_length / actual_length;

            //

            if (std::get<1>(it->second)->isFixed() == true && std::get<0>(it->second)->isFixed() == false) {
                glm::vec3 actual_length_vector = actual_length * glm::normalize(
                        std::get<0>(it->second)->getPosition() - std::get<1>(it->second)->getPosition());
                std::get<0>(it->second)->setPosition(actual_length_vector * correction_factor_x +
                                                     (std::get<1>(it->second)->getPosition()));
            }
            if (std::get<0>(it->second)->isFixed() == true && std::get<1>(it->second)->isFixed() == false) {
                glm::vec3 actual_length_vector = actual_length * glm::normalize(
                        std::get<1>(it->second)->getPosition() - std::get<0>(it->second)->getPosition());
                std::get<1>(it->second)->setPosition(actual_length_vector * correction_factor_x +
                                                     (std::get<0>(it->second)->getPosition()));
            }

            if (std::get<0>(it->second)->isFixed() == false && std::get<1>(it->second)->isFixed() == false) {
                glm::vec3 actual_length_vector = actual_length * glm::normalize(
                        std::get<1>(it->second)->getPosition() - std::get<0>(it->second)->getPosition());
                glm::vec3 shortened_vector = (actual_length_vector - actual_length_vector * correction_factor_x) / 2.0f;
                std::get<0>(it->second)->setPosition(shortened_vector + (std::get<0>(it->second)->getPosition()));
                std::get<1>(it->second)->setPosition(std::get<1>(it->second)->getPosition() - shortened_vector);

            }
            it++;
        }
    }
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            //Update original position and volocity
            getParticle(i,j).setOriginalPosition(getParticle(i,j).getPosition());
            getParticle(i,j).setOriginalVelocity(getParticle(i,j).getVelocity());
        }
    }




    if(initial == true){
        initial =false;
    }
    // commented out because an animated bounding box can be weird
    // computeBoundingBox();
}

std::vector<std::tuple<ClothParticle *, float, int>>
Cloth::identify_strings(int i, int j, std::vector<std::pair<ClothParticle *, float>> &structural,
                        std::vector<std::pair<ClothParticle *, float>> &shear,
                        std::vector<std::pair<ClothParticle *, float>> &flexion,
                        std::unordered_map< std::string , std::tuple<ClothParticle*, ClothParticle*,int, float > > &springs) {
    collect_springs(i, j, structural, shear, flexion);

    //combine all strings
    std::vector<std::tuple<ClothParticle*, float, int>> all_springs;
    for (int l = 0; l < structural.size(); ++l) {
        all_springs.push_back(std::make_tuple(structural[l].first, structural[l].second, 0));
        if(initial == true){
            glm::vec3 start = structural[l].first->getPosition();
            glm::vec3 end = getParticle(i,j).getPosition();
            if(springs.find(spring_to_string(start,end)) == springs.end() &&
               springs.find(spring_to_string(end,start)) == springs.end())
            {
                std::get<0>(springs[spring_to_string(start,end)]) = structural[l].first;
                std::get<1>(springs[spring_to_string(start,end)]) = &getParticle(i,j);
                std::get<2>(springs[spring_to_string(start,end)]) = 0;
                std::get<3>(springs[spring_to_string(start,end)]) = structural[l].second;
            }
        }
    }
    for (int l = 0; l < shear.size(); ++l) {
        all_springs.push_back(std::make_tuple(shear[l].first, shear[l].second, 1));
        if(initial == true) {
            glm::vec3 start = shear[l].first->getPosition();
            glm::vec3 end = getParticle(i, j).getPosition();
            if (springs.find(spring_to_string(start, end)) == springs.end() &&
                springs.find(spring_to_string(end, start)) == springs.end()) {
                std::get<0>(springs[spring_to_string(start, end)]) = shear[l].first;
                std::get<1>(springs[spring_to_string(start, end)]) = &getParticle(i, j);
                std::get<2>(springs[spring_to_string(start, end)]) = 1;
                std::get<3>(springs[spring_to_string(start, end)]) = shear[l].second;
            }
        }
    }
    for (int l = 0; l < flexion.size(); ++l) {
        all_springs.push_back(std::make_tuple(flexion[l].first, flexion[l].second, 2));
    }
    return all_springs;
}

std::string Cloth::spring_to_string(glm::vec3 a, glm::vec3 b){
    std::string s;
    s = "("+std::to_string(a[0])+","+std::to_string(a[1])+","+std::to_string(a[2])+")" +
            "("+std::to_string(b[0])+","+std::to_string(b[1])+","+std::to_string(b[2])+")";
    return s;
}

void Cloth::collect_springs(int i, int j, std::vector<std::pair<ClothParticle *, float>> &structural,
                            std::vector<std::pair<ClothParticle *, float>> &shear,
                            std::vector<std::pair<ClothParticle *, float>> &flexion) {//Collect Structural
    //Collect Structual
    if(i+1 < nx)
    {
        structural.push_back( std::make_pair(&getParticle(i + 1, j), natural_x) );
    }
    if(i-1 >= 0)
    {
        structural.push_back( std::make_pair(&getParticle(i - 1, j), natural_x) );
    }
    if(j+1 < ny)
    {
        structural.push_back( std::make_pair(&getParticle(i , j + 1), natural_y) );
    }
    if(j-1 >= 0)
    {
        structural.push_back( std::make_pair(&getParticle(i , j - 1), natural_y) );
    }


    //Collect Shear
    if((i + 1 < nx) && (j + 1 < ny))
    {
        shear.push_back( std::make_pair(&getParticle(i + 1 , j + 1), natural_xy));
    }
    if((i + 1 < nx) && (j - 1 >= 0))
    {
        shear.push_back( std::make_pair(&getParticle(i + 1 , j - 1), natural_xy));
    }
    if((i - 1 >= 0) && (j + 1 < ny))
    {
        shear.push_back( std::make_pair(&getParticle(i - 1 , j + 1), natural_xy));
    }
    if((i - 1 >= 0) && (j - 1 >= 0))
    {
        shear.push_back( std::make_pair(&getParticle(i - 1 , j - 1), natural_xy));
    }


    //Collect flexion
    if(i + 2 < nx)
    {
        flexion.push_back( std::make_pair(&getParticle(i + 2 , j), 2.0f * natural_x));
    }
    if(i - 2 >= 0)
    {
        flexion.push_back(std::make_pair(&getParticle(i - 2 , j), 2.0f * natural_x));
    }
    if(j - 2 >= 0)
    {
        flexion.push_back(std::make_pair(&getParticle(i , j - 2), 2.0f * natural_y));
    }
    if(j + 2 < ny)
    {
        flexion.push_back(std::make_pair(&getParticle(i , j + 2), 2.0f * natural_y));
    }
}

