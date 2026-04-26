#include "rsolverwave.h"

RSolverWave::RSolverWave(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData)
    : RSolverGeneric(pModel,modelFileName,convergenceFileName,sharedData)
{
    this->problemType = R_PROBLEM_WAVE;
}

RSolverWave::~RSolverWave()
{

}

bool RSolverWave::hasConverged() const
{
    return true;
}

void RSolverWave::updateScales()
{

}

void RSolverWave::recover()
{

}

void RSolverWave::prepare()
{
    RBVector waveDisplacementExplicitFlags;

    this->generateNodeBook(R_PROBLEM_WAVE);
    this->generateVariableVector(R_VARIABLE_WAVE_DISPLACEMENT,this->elementWaveDisplacement,waveDisplacementExplicitFlags,true,this->firstRun,this->firstRun);
//    this->generateVariableVector(R_VARIABLE_W,this->elementWaveDisplacement,waveDisplacementExplicitFlags,true,firstTime,firstTime);

    this->b.resize(this->nodeBook.getNEnabled());
    this->x.resize(this->nodeBook.getNEnabled());

    this->A.clear();
    this->b.fill(0.0);
    this->x.fill(0.0);
}

void RSolverWave::solve()
{

}

void RSolverWave::process()
{

}

void RSolverWave::store()
{

}

void RSolverWave::statistics()
{

}
