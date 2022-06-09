/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::Testing
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// This module defines the %Testing class that computes for each vertex of a
/// triangulation the average scalar value of itself and its direct neighbors.
///
/// \b Related \b publication: \n
/// 'Testing'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2021.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

namespace ttk {

  /**
   * The Testing class provides methods to compute for each vertex of a
   * triangulation the average scalar value of itself and its direct neighbors.
   */
  class Testing : virtual public Debug {

  public:
    Testing();

    inline void
      setOutput(std::vector<std::pair<ttk::SimplexId, char>> *criticalPoints) {
      criticalPoints_ = criticalPoints;
    }

    /**
     * TODO 2: This method preconditions the triangulation for all operations
     *         the algorithm of this module requires. For instance,
     *         preconditionVertexNeighbors, preconditionBoundaryEdges, ...
     *
     *         Note: If the algorithm does not require a triangulation then
     *               this method can be deleted.
     */
    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      return triangulation->preconditionVertexNeighbors();
    }

    /**
     * TODO 3: Implmentation of the algorithm.
     *
     *         Note: If the algorithm requires a triangulation then this
     *               method must be called after the triangulation has been
     *               preconditioned for the upcoming operations.
     */
    template <class dataType,
              class triangulationType = ttk::AbstractTriangulation>
    int computeCriticalPoints(const dataType *inputData,
                        const triangulationType *triangulation) const {
      // start global timer
      ttk::Timer globalTimer;

      // print horizontal separator
      this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator

      // print input parameters in table format
      this->printMsg({
        {"#Threads", std::to_string(this->threadNumber_)},
        {"#Vertices", std::to_string(triangulation->getNumberOfVertices())},
      });
      this->printMsg(ttk::debug::Separator::L1);

      // -----------------------------------------------------------------------
      // Compute Vertex Averages
      // -----------------------------------------------------------------------
      {

        // compute the average of each vertex in parallel
        ttk::SimplexId nVertices = static_cast<ttk::SimplexId>(triangulation->getNumberOfVertices());

        this->printMsg("I was at least here");

// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(this->threadNumber_)
// #endif
        for(ttk::SimplexId i = 0; i < nVertices; i++) {
          double vi = inputData[3*i];
          double vj = inputData[3*i + 1];
          if (isZero(vi) && isZero(vj)) {
            criticalPoints_->push_back({i, (char) 3});
          }
          else {
            // Check if critical point is inside a cell
            size_t nNeighbors = triangulation->getVertexNeighborNumber(i);
            ttk::SimplexId neighborId;
            int isCriticalPoint = (isZero(vi, 1e-1) && isZero(vj, 1e-1));
            for(size_t j = 0; j < nNeighbors; j++) {
              if (isCriticalPoint == 0) break;
              
              triangulation->getVertexNeighbor(i, j, neighborId);
              double ni = inputData[3 * neighborId];
              double nj = inputData[3 * neighborId + 1];
              if (!(vi < ni && vj < nj)) {
                isCriticalPoint = 0;
              }
            }
            if (isCriticalPoint) {
              criticalPoints_->push_back({i, (char) 3});
            }
          }
        }

        this->printMsg("I was here too");
      }

      this->printMsg("Not the end");

      // ---------------------------------------------------------------------
      // print global performance
      // ---------------------------------------------------------------------
      {
        this->printMsg(ttk::debug::Separator::L2); // horizontal '-' separator
        this->printMsg(
          "Complete", 1, globalTimer.getElapsedTime() // global progress, time
        );
        this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
      }

      this->printMsg("End");


      return 1; // return success
    }
  protected:
    std::vector<std::pair<ttk::SimplexId, char>> * criticalPoints_;

  private:
    int isZero(const double value, const double thresh=1e-4) const {
      return (value <= thresh && value >= -thresh) ? 1 : 0;
    }
    

  }; // Testing class

} // namespace ttk
