/*
 * NOTE:
 *      对于粘性热耗散项，需要考虑壁面附近粘性加热的影响，这与壁面函数相关。对于
 *  两方程模型(如k-epsilon/omega)，基于Jayatilleke模型并对OpenFOAM稍作修改基本可
 *  以满足要求；但是对于一方程的SA模型，OpenFOAM目前的算法不能满足要求。
 *      需要改进的方向：1，对OpenFOAM的alphat壁面函数进行重构，在更高的层次提供
 *  计算热耗散通量的API；2，为SA模型提供合适的alphat壁面函数。
 */
surfaceScalarField D
(
    "D",
    fvc::dotInterpolate(mesh.Sf(), turbulence->devRhoReff()) & fvc::interpolate(U)
);

const volScalarField& alphat = turbulence->alphat();
surfaceScalarField::Boundary& Dbnd = D.boundaryFieldRef();

const fvBoundaryMesh& meshBnd = mesh.boundary();

forAll(meshBnd, patchi)
{
    const fvPatch& patch = meshBnd[patchi];
    using namespace compressible;

    const fvPatchScalarField& alphatp = alphat.boundaryField()[patchi];

    if
    (
        isA<wallFvPatch>(patch)
     && isA<alphatWallFunctionFvPatchScalarField>(alphatp)
    )
    {
        Dbnd[patchi] =
          - dynamic_cast<const alphatWallFunctionFvPatchScalarField&>
            (
                alphatp
            ).qD()*patch.magSf();
    }
}
