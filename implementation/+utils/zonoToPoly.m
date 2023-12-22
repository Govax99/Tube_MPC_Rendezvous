function P = zonoToPoly(obj)
    % Plot a constrained zonotope
    % Optional arguments:
    % dims: vector with dimentions to project and plot the object
    % color: color definition of the interior of the object to be plotted
    % alpha: transparency of the color defined
    % Default color: red. Default alpha: 1. Default dims: [1 2 3], [1 2], [1]
    % for 3-, 2- and 1-dimensional objects respectively.
    %
    % Usage options:
    % obj.plot([1 2])
    % obj.plot('b',0.1)
    % obj.plot([1 2],'b',0.1);
    %
    c = obj.c;
    G = obj.G;
    
    Box = Polyhedron('lb',-ones(obj.nG,1),'ub',ones(obj.nG,1));
    P = c + G*Box;
end

