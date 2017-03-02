
double Fluid::AdjustForIncompressibility() {


    //max divergence
    double max_div = 0;


    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                Cell *c = getCell(i,j,k);
                if (c->numParticles()!=0) {
                    //defalut divider is 6, is any surface is unavaliable, decrease divider
                    double divider = 6;
                    //calculate divergence
                    double divergence = (get_new_u_plus(i, j, k) - get_new_u_plus(i - 1, j, k)) +
                                        (get_new_v_plus(i, j, k) - get_new_v_plus(i, j - 1, k)) +
                                        (get_new_w_plus(i, j, k) - get_new_w_plus(i, j, k - 1));

                    if (i == nx - 1) {
                        divider-=1;
                    };
                    if (j == ny - 1) {
                        divider-=1;
                    };
                    if (k == nz - 1){
                        divider-=1;
                    };
                    if (i == 0) {
                        divider-=1;
                    };
                    if (j == 0) {
                        divider-=1;
                    };
                    if (k == 0) {
                        divider-=1;
                    };

                    //the volocity portion to adjust in every direction
                    double div_delta = divergence / divider;

                    if (i != nx - 1) {
                        set_new_u_plus(i, j, k, get_new_u_plus(i, j, k) - div_delta);
                    };
                    if (j != ny - 1) {
                        set_new_v_plus(i, j, k, get_new_v_plus(i, j, k) - div_delta);
                    };
                    if (k != nz - 1){
                        set_new_w_plus(i, j, k, get_new_w_plus(i, j, k) - div_delta);
                    };
                    if (i != 0) {
                        set_new_u_plus(i - 1, j, k, get_new_u_plus(i - 1, j, k) + div_delta);
                    };
                    if (j != 0) {
                        set_new_v_plus(i, j - 1, k, get_new_v_plus(i, j - 1, k) + div_delta);
                    };
                    if (k != 0) {
                        set_new_w_plus(i, j, k - 1, get_new_w_plus(i, j, k - 1) + div_delta);
                    };
                }
            }
        }
    }


    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                Cell *c = getCell(i,j,k);
                if (c->numParticles()!=0) {
                    double divergence =
                            fabs((1 / dx) * (get_new_u_plus(i, j, k) - get_new_u_plus(i - 1, j, k)) +
                                 (1 / dy) * (get_new_v_plus(i, j, k) - get_new_v_plus(i, j - 1, k)) +
                                 (1 / dz) * (get_new_w_plus(i, j, k) - get_new_w_plus(i, j, k - 1)));
                    if (divergence > max_div) {
                        max_div = divergence;
                    }
                }
            }
        }
    }
    // return the maximum divergence
    // (will iterate for specified # of iterations or until divergence is near zero)
    return max_div;
}

// ==============================================================





glm::vec3 Fluid::getInterpolatedVelocity(const glm::vec3 &pos) const {

    int i = int(floor(pos.x/dx)); if (i < 0) i = 0; if (i >= nx) i = nx-1;
    int j = int(floor(pos.y/dy)); if (j < 0) j = 0; if (j >= ny) j = ny-1;
    int k = int(floor(pos.z/dz)); if (k < 0) k = 0; if (k >= nz) k = nz-1;
    glm::vec3 cube_center = glm::vec3((0.5 + float(i)) * dx, (0.5 + float(j)) * dy, (0.5 + float(k)) * dz);
    float cell_volume = dx * dy * dz;
    float u_component = 0;
    float v_component = 0;
    float w_component = 0;


    for (int t =i-1 ; t<=i+1 ; ++t){
        for (int p = j-1; p <= j+1; ++p) {
            for (int s = k-1; s <= k+1 ; ++s) {
                //deal with u direction
                glm::vec3 u_current_cube = glm::vec3((1.0f + float(t)) * dx, (0.5f + float(p)) * dy,  (0.5f + float(s)) * dz);
                if((fabs(u_current_cube[0] - pos.x) < dx) &&
                   (fabs(u_current_cube[1] - pos.y) < dy) &&
                   (fabs(u_current_cube[2] - pos.z) < dz)){
                    //if one of the u we want
                    float u_speed = get_u_plus(t,　p,s);
                    float overlapVolume = (dx - fabs(pos.x - u_current_cube[0])) *
                                          (dy - fabs(pos.y - u_current_cube[1])) *
                                          (dz - fabs(pos.z - u_current_cube[2]));
                    u_component += u_speed * (overlapVolume / cell_volume);
                }

                //deal with v direction
                glm::vec3 v_current_cube = glm::vec3((0.5f + float(t)) * dx, (1.0f + float(p)) * dy,  (0.5f + float(s)) * dz);
                if((fabs(v_current_cube[0] - pos.x) < dx) &&
                   (fabs(v_current_cube[1] - pos.y) < dy) &&
                   (fabs(v_current_cube[2] - pos.z) < dz)){
                    //if one of the v we want
                    float v_speed = get_v_plus(t,p,s);
                    float overlapVolume = (dx - fabs(pos.x - v_current_cube[0])) *
                                          (dy - fabs(pos.y - v_current_cube[1])) *
                                          (dz - fabs(pos.z - v_current_cube[2]));
                    v_component += v_speed * (overlapVolume / cell_volume);
                }

                //deal with w direction
                glm::vec3 w_current_cube = glm::vec3((0.5f + float(t)) * dx, (0.5f + float(p)) * dy,  (1.0f + float(s)) * dz);
                if((fabs(w_current_cube[0] - pos.x) < dx) &&
                   (fabs(w_current_cube[1] - pos.y) < dy) &&
                   (fabs(w_current_cube[2] - pos.z) < dz)){
                    //if one of the w we want
                    float w_speed = get_w_plus(t,p,s);
                    float overlapVolume = (dx -fabs(pos.x - w_current_cube[0])) *
                                          (dy - fabs(pos.y - w_current_cube[1])) *
                                          (dz - fabs(pos.z - w_current_cube[2]));
                    w_component += w_speed * (overlapVolume / cell_volume);
                }

            }
        }
    }

    return glm::vec3(u_component,v_component,w_component);
    //
    // *********************************************************************


}


float Fluid::calculate_cell_volume(glm::vec3 front_bottomn_left, glm::vec3 back_top_right) {
    float volume;
    volume = fabsf(front_bottomn_left[0] - back_top_right[0]) *
             fabsf(front_bottomn_left[1] - back_top_right[1]) *
             fabsf(front_bottomn_left[2] - back_top_right[2]);
    return volume;

}
