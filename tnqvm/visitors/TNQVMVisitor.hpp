/***********************************************************************************
 * Copyright (c) 2018, UT-Battelle
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *   * Neither the name of the xacc nor the
 *     names of its contributors may be used to endorse or promote products
 *     derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Contributors:
 *   Initial API and implementation - Alex McCaskey
 *
 **********************************************************************************/
#ifndef TNQVM_TNQVMVISITOR_HPP_
#define TNQVM_TNQVMVISITOR_HPP_

#include "Identifiable.hpp"
#include "AllGateVisitor.hpp"
#include "xacc.hpp"

using namespace xacc;
using namespace xacc::quantum;

namespace tnqvm {
// Struct described an observable term (Pauli operators) for expectation value calculation.
// i.e. get the result of `<psi| Hamiltonian |psi>` for an arbitrary Hamiltonian.
// The Hamiltonian must be expressed as a sum of products of gates (with a complex coefficient for each term):
// e.g. a1*Z1Z2 + a2*X0X1 + .. + Y2Y3, etc.
struct ObservableTerm
{
    ObservableTerm(const std::vector<std::shared_ptr<Instruction>>& in_operatorsInProduct, const std::complex<double>& in_coeff = 1.0):
      coefficient(in_coeff),
      operators(in_operatorsInProduct)
    {}

    std::complex<double> coefficient;
    std::vector<std::shared_ptr<Instruction>> operators;
};
// Typedef for the whole Hamiltonian, i.e. a vector of ObservableTerms
using ObservableExpr =  std::vector<ObservableTerm>;    

class TNQVMVisitor : public AllGateVisitor, public OptionsProvider,
                     public xacc::Cloneable<TNQVMVisitor> {
public:
  virtual void initialize(std::shared_ptr<AcceleratorBuffer> buffer, int nbShots = 1) = 0;
  virtual const double
  getExpectationValueZ(std::shared_ptr<CompositeInstruction> function) = 0;

  // Compute the expectation value of an observable (e.g. energy from Hamiltonian)
  // Tensor processing backends of TNQVM can use tensor contraction to compute this value efficiently (in terms of memory)
  // Note: not all TNQVM backend has method implemented.
  virtual double getExpectationValue(std::shared_ptr<AcceleratorBuffer>& buffer, std::shared_ptr<CompositeInstruction>& function, const ObservableExpr& in_observable) {
    throw std::runtime_error(this->name() + " backend does not support direct observable expectation calculation.\n");
  }

  virtual const std::vector<std::complex<double>> getState() {
    return std::vector<std::complex<double>>{};
  }
  virtual void finalize() = 0;
 
protected:
  std::shared_ptr<AcceleratorBuffer> buffer;
};

} // namespace tnqvm

#endif /* TNQVM_TNQVMVISITOR_HPP_ */
