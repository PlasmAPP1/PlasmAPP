```MATLAB
function interp = InterpolationF(position,N_particles,aux_vector)
    switch Interpolation_method
        case 'Nearest Grid Point (NGP)'
                if N_particles == 0
                    interp = 0;
                else
                    % Project Particles to grid (Subscripts denote the beams):
                    im = floor(position/dx)+1;
                    ip = im + 1;
                    % im: First node. ip = i+1.
                    Project = [im;ip];                           % Concatenate nodes which particle is between.

                    % Boundary conditions on the projection:
                    Project = PeriodicBC(,Project,1,Ng);

                    % Compute what fraction of the particle size that lies on the two nearest cells:
                    Fim = abs((position/dx) - im);       % Fim: Fraction of the particle that lies in the node i.
                    Fip = 1 - Fim;                                   % Fip: Fraction of the particle that lies in the node i1.
                    Fraction = [Fim;Fip];                            % Concatenate Fim and Fip.
                    Fraction(Fraction>0.5) = 1;
                    Fraction(Fraction<0.5) = 0;

                    % Interpolation matrix: Gives me a matrix of N1*Ng, were it writes the fraction of the particles in the position
                    % (p1,Project_1).
                    interp = sparse(aux_vector, Project, Fraction, N_particles, Ng);   % Matrix of N1*Ng, p1 : position i, Project_1 : position j.
                end
        case 'Cloud in Cell (CIC)'
                if N_particles == 0
                    interp = 0;
                else
                    % Project Particles to grid (Subscripts denote the beams):
                    im = floor(position/dx-.5) + 1;
                    ip = im + 1;
                    % im: First node. ip = i+1.
                    Project = [im;ip];                           % Concatenate nodes which particle is between.

                    % Boundary conditions on the projection:
                    Project = PeriodicBC(,Project,1,Ng);

                    % Compute what fraction of the particle size that lies on the two nearest cells:
                    Fim = abs((position/dx)- im+.5);       % Fim: Fraction of the particle that lies in the node i.
                    Fip = 1 - Fim;                                   % Fip: Fraction of the particle that lies in the node i1.
                    Fraction = [Fim;Fip];                            % Concatenate Fim and Fip.

                    % Interpolation matrix: Gives me a matrix of N1*Ng, were it writes the fraction of the particles in the position
                    % (p1,Project_1).
                    interp = sparse(aux_vector, Project, Fraction, N_particles, Ng);   % Matrix of N1*Ng, p1 : position i, Project_1 : position j.
                end
    end         
end