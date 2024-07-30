% dv e dw tem tamanho N2x1
function prediction = modelo(x_k, y_k, theta_k, v_kmin1, w_kmin1, dv, dw, eta_x_k, eta_y_k, eta_theta_k, N2, Ts)
    x(1) = x_k;
    y(1) = y_k;
    theta(1) = theta_k;

    v = v_kmin1*ones(N2,1) + tril(ones(N2))*dv;
    w = w_kmin1*ones(N2,1) + tril(ones(N2))*dw;

    for k=1:N2
        x(k+1) = x(k) + Ts*v(k)*cos(theta(k));
        y(k+1) = y(k) + Ts*v(k)*sin(theta(k));
        theta(k+1) = theta(k) + Ts*w(k);
    end

    prediction(:,1) = x(2:end)' + eta_x_k;
    prediction(:,2) = y(2:end)' + eta_y_k;
    prediction(:,3) = theta(2:end)' + eta_theta_k;
end