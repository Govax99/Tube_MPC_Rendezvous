function plotSolutionTrajectoryGeneration(trajectoryGenerator)
    % DESCRIPTION: Plot trajectory given the results from the generator.
    % 1) propagate tumbler obtain end attitude
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [~, s_tumbler] = ode45(@(t,s) dynamics.odetumbler(t,s,trajectoryGenerator.context.parameters), ...
        [0, trajectoryGenerator.context.ref.maneuverTime], trajectoryGenerator.context.ref.initialTumblerState, options);
    s_tumbler = s_tumbler';
    config.q = s_tumbler(4:7,:);
    wT_fin = s_tumbler(1:3,end);
    qT_fin = s_tumbler(4:7,end);
    
    % 2) obtain final position and velocity for rendezvous
    pBerth = trajectoryGenerator.context.parameters.pBerth;
    R_LT = quat.quat2rotm(qT_fin);
    p_final = R_LT*pBerth;
    w_final_L = R_LT*wT_fin - [0; 0; trajectoryGenerator.context.parameters.OM];
    v_final = cross(w_final_L, p_final);
    
    t_guess = trajectoryGenerator.problem.guess.time;
    p_guess = trajectoryGenerator.problem.guess.state(1:3,:);
    v_guess = trajectoryGenerator.problem.guess.state(4:6,:);
    
    t = trajectoryGenerator.solutionInertialFrame.t;
    q = trajectoryGenerator.solutionInertialFrame.z(1:3,:);
    dq = trajectoryGenerator.solutionInertialFrame.z(4:6,:);
    %u = soln(end).grid.control;
    
    sat = visualization.Satellite(config);
    title("Satellite Rendezvous - LVLH frame")
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('z [m]')
    axis('equal')
    lim = 15;
    xlim([-lim, lim])
    ylim([-lim, lim])
    zlim([-lim, lim])
    grid on;
    hold on;
    view(45,45)
    nTimeGuess = length(s_tumbler);
    sat.plot(nTimeGuess)
    plot3([0 p_final(1)], [0 p_final(2)], [0 p_final(3)])
    plot3([0 w_final_L(1)], [0 w_final_L(2)], [0 w_final_L(3)])
    plot3([p_final(1) p_final(1)+v_final(1)], [p_final(2) p_final(2)+v_final(2)], [p_final(3) p_final(3)+v_final(3)])
    
    
    guess = quiver3(p_guess(1,:),p_guess(2,:),p_guess(3,:),v_guess(1,:),v_guess(2,:),v_guess(3,:));
    hold on;
    plot3(p_guess(1,:),p_guess(2,:),p_guess(3,:))
    
    sol = quiver3(q(1,:),q(2,:),q(3,:),dq(1,:),dq(2,:),dq(3,:));
    hold on;
    plot3(q(1,:),q(2,:),q(3,:),'--')
%     legend([guess, sol],["Guess trajectory","Solution trajectory"])


    figure
    plot(t,q)
    title("Positions in Inertial frame")
    xlabel('t [s]')
    ylabel('Positions [m]')

    figure
    plot(t,dq)
    title("Velocities in Inertial frame")
    xlabel('t [s]')
    ylabel('Velocities [m/s]')
end
