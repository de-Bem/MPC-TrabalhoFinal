function [Hw,fw,Aw,b] = gerenciador_restricao(G_total,resposta_livre,N1,N2,Nuv,Nuw,n_epsilons,v_kmin1,w_kmin1,aux_Hw,aux_fw,aux_Aw,A_restr_controle,A_restr_increm,b_aux_restr_saida,b_aux_restr_controle,b_restr_increm,Hqp,fqp)    

    A_restr_saida = [ G_total; -G_total];

    b_restr_saida = b_aux_restr_saida + ...
        [-resposta_livre(N1:N2,1);-resposta_livre(N1:N2,2);-resposta_livre(N1:N2,3);
          resposta_livre(N1:N2,1); resposta_livre(N1:N2,2); resposta_livre(N1:N2,3)];

    b_restr_controle = b_aux_restr_controle + ...
        [-v_kmin1*ones(Nuv,1); v_kmin1*ones(Nuv,1);
         -w_kmin1*ones(Nuw,1); w_kmin1*ones(Nuw,1)];

    A = [A_restr_saida;
         A_restr_controle;
         A_restr_increm];

    b = [b_restr_saida;
         b_restr_controle;
         b_restr_increm];

    % variável de folga:

    aux_Hw(1:end-n_epsilons,1:end-n_epsilons) = Hqp; Hw=aux_Hw;
    aux_fw(1:end-n_epsilons) = fqp; fw=aux_fw;
    aux_Aw(:,1:end-n_epsilons) = A; Aw=aux_Aw;

end