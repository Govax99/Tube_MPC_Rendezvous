% launch with: runtests('rotateTest')
function tests = rotateTest
    tests = functiontests(localfunctions);
end

function setup(testCase)
    % Add necessary setup code here
end

function setupOnce(testCase)  % do not change function name
% set a new path, for example
    cd("..\..\")
end

function teardown(testCase)
    % Add necessary cleanup code here
end

function testRotateIdentityQuaternion(testCase)
    % Test that rotating a vector with an identity quaternion results in no change
    v = [1; 2; 3];
    q = [0; 0; 0; 1];  % Identity quaternion
    vp = quat.rotate(v, q);
    expected_vp = v;
    assert(isequal(vp, expected_vp));
end

function testRotate90DegreesAboutXAxis(testCase)
    % Test rotating a vector 90 degrees about the X-axis
    v = [1; 0; 0];
    q = [1; 0; 0; 0];  % 90-degree rotation about X-axis
    vp = quat.rotate(v, q);
    expected_vp = [1; 0; 0];
    assert(isequaln(vp, expected_vp));
end

function testRotate45DegreesAboutZAxis(testCase)
    % Test rotating a vector 90 degrees about the X-axis
    v = [1; 0; 0];
    q = [ 0, 0, 0.3826834, 0.9238795 ];  % 90-degree rotation about X-axis
    vp = quat.rotate(v, q);
    expected_vp = [cosd(45); sind(45); 0];
    assert(utils.isalmost(vp, expected_vp, 1e-6));
end

function testRotateArbitraryVector(testCase)
    % Test rotating an arbitrary vector with a non-identity quaternion
    N = 1000;
    for i = 1:N
        qMat = compact(randrot);
        q = [qMat(2:4), qMat(1)];
        v = rand([1,3]);
        vp = quat.rotate(v, quat.conj(q));
        % There are 2 main differences between my quaternions and matlab
        % quaternions: 1) [qv, qs] vs [qs, qv] 2) mine is active rotations,
        % theirs is passive rotations
        vpMat = quaternionRotationMatlab(v, qMat);
        assert(utils.isalmost(vp, vpMat));
    end
end

function testRotateAgainstMatrixProduct(testCase)
    N = 1000;
    for i = 1:N
        qMat = compact(randrot);
        q = [qMat(2:4), qMat(1)];
        R = quat.quat2rotm(q);
        v = rand([1,3]);
        vp = quat.rotate(v, q);
        % There are 2 main differences between my quaternions and matlab
        % quaternions: 1) [qv, qs] vs [qs, qv] 2) mine is active rotations,
        % theirs is passive rotations
        vpMat = R*v';
        assert(utils.isalmost(vp, vpMat));
    end
end

function vp = quaternionRotationMatlab(v, q)
    vp = quatrotate(q, v);
    vp = vp';
end

