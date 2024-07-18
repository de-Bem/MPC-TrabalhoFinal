function [Ref,Qe] = gerenciador_referencia(N,ref_idx,delta_x,delta_y,delta_theta,Referencia_x,Referencia_y,Referencia_theta)

    Ref_x = Referencia_x(ref_idx)*ones(N,1);
    Ref_y = Referencia_y(ref_idx)*ones(N,1);
    Ref_theta = Referencia_theta(ref_idx)*ones(N,1);
    Ref = [Ref_x; Ref_y; Ref_theta];
    
    dx = Referencia_x(ref_idx)-Referencia_x(ref_idx-1);
    dy = Referencia_y(ref_idx)-Referencia_y(ref_idx-1);
    dtheta = Referencia_theta(ref_idx)-Referencia_theta(ref_idx-1);

    if dx == 0
        dx = 1;
    end
    if dy == 0
        dy = 1;
    end
    if dtheta == 0
        dtheta = 1;
    end

    Qe_x = delta_x*eye(N)/(dx^2*N); 
    Qe_y = delta_y*eye(N)/(dy^2*N);
    Qe_theta = delta_theta*eye(N)/(dtheta^2*N);
    Qe = [    Qe_x   zeros(N)  zeros(N);
            zeros(N)   Qe_y    zeros(N);
            zeros(N) zeros(N) Qe_theta];

end