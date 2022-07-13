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
      int result = triangulation->preconditionCellTriangles();
      return result;
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

      // // print input parameters in table format
      // this->printMsg({
      //   {"#Threads", std::to_string(this->threadNumber_)},
      //   {"#Vertices", std::to_string(triangulation->getNumberOfVertices())},
      // });
      // this->printMsg(ttk::debug::Separator::L1);

      // -----------------------------------------------------------------------
      // Compute Vertex Averages
      // -----------------------------------------------------------------------
      
      // criticalPoints_->push_back({0, (char) 0});

      this->printMsg("Beginning cell checking");
      this->printMsg("Number of cells: " +  std::to_string(triangulation->getNumberOfCells()));
      SimplexId nCells = triangulation->getNumberOfCells();

      criticalPoints_->clear();

      // int n0 = 0, n1 = 0;
      
      // for (SimplexId i = 40048; i < 40050; i++) {
      for (SimplexId i = 0; i < nCells; i++) {
        SimplexId nVertices = triangulation->getCellVertexNumber(i);
        // this->printMsg("Vertices in Cell #" + std::to_string(i) + ": " + std::to_string(nVertices));
        // this->printMsg(" Cell #" + std::to_string(i));
        assert(nVertices == 3);
        double values[nVertices * 2];
        SimplexId vertices[nVertices];

        for (SimplexId j = 0; j < nVertices; ++j) {
          SimplexId vertexId;
          triangulation->getCellVertex(i, j, vertexId);

          // this->printMsg("Vertex #" + std::to_string(vertexId));
          values[j * 2] = inputData[3 * vertexId];
          values[j * 2 + 1] = inputData[3 * vertexId + 1];
          vertices[j] = vertexId;
          // this->printMsg("Vertex: " + std::to_string(j) +  ", " +  std::to_string(vertexId) +  ", " +  
          //   "{" + std::to_string(vi) + ", " + std::to_string(vj) + "}");
        }

        if (!zeroInBoundingBox(values[0], values[1], values[2], values[3], values[4], values[5])) {
          continue;
        }

        int ret = zeroInTriangle(values[0], values[1], values[2], values[3], values[4], values[5]);
        if (ret >= 0) {
          criticalPoints_->push_back({vertices[0], (char) ret});
          // this->printMsg("Simplex Id: " + std::to_string(i));
        }
      }

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

    double det(const double vi1, const double vj1, const double vi2, const double vj2, const double vi3, const double vj3) const {
      return vi1 * (vj2 - vj3) + (vi2 * vj3 - vj2 * vi3) - vj1 * (vi2 - vi3);
    }

    inline int sign(const double x) const {
      if (x < 0) return -1;
      if (x > 0) return 1;
      return 0;
    }
    
    int positive(const double vi1, const double vj1, const double vi2, const double vj2, const double vi3, const double vj3) const {
      double val = det(vi1, vj1, vi2, vj2, vi3, vj3);
      if (val > 0) return 1;
      if (val < 0) return -1;
      // this->printMsg("Zero determinant!!");
      // this->printMsg("Values: "   
      //       "{" + std::to_string(vi1) + ", " + std::to_string(vj1) + "}" + 
      //       "{" + std::to_string(vi2) + ", " + std::to_string(vj2) + "}" + 
      //       "{" + std::to_string(vi3) + ", " + std::to_string(vj3) + "}" 
      //       );
      // -P(2, 1)                      E(1, 2) 
      if (vi2 != 0.0) return -sign(vi2);
      // +P(3, 1)                      E(1, 2) 
      if (vi3 != 0.0) return sign(vi3);
      // -P(3, 2)                      E(1, 1) 
      if (vj3 != 0.0) return -sign(vj3);
      // +P(2, 2)                      E(1, 1) 
      if (vj2 != 0.0) return sign(vj2);
      // -P(3, 1)                      E(2, 2) 
      if (vi3 != 0.0) return -sign(vi3);
      // +P(1, 1)                      E(2, 2) 
      if (vi1 != 0.0) return sign(vi1);
      // +                             E(1, 1) E(2, 2) 
      return 1;
    }

    int zeroInTriangle(const double vi1, const double vj1, const double vi2, const double vj2, const double vi3, const double vj3) const {
      
      int sign = positive(vi1, vj1, vi2, vj2, vi3, vj3);
      // if (sign == 0) return 0;
      int sign1 = positive(vi1, vj1, vi2, vj2, 0, 0);
      if (sign1 != sign) return -1;
      // if (sign1 == 0) return 0;
      int sign2 = -positive(vi1, vj1, vi3, vj3, 0, 0);
      if (sign2 != sign) return -1;
      // if (sign2 == 0) return 0;
      int sign3 = positive(vi2, vj2, vi3, vj3, 0, 0);
      if (sign3 != sign) return -1;
      return 1;
    }

    inline bool zeroInBoundingBox(const double vi1, const double vj1, const double vi2, const double vj2, const double vi3, const double vj3) const {
      double maxX = std::max(vi1, std::max(vi2, vi3));
      double maxY = std::max(vj1, std::max(vj2, vj3));
      double minX = std::min(vi1, std::min(vi2, vi3));
      double minY = std::min(vj1, std::min(vj2, vj3));
      if (maxX >= 0.0 && maxY >= 0.0 && minY <= 0.0 && minX <= 0.0) {
        return true;
      }
      return false;
    }

  }; // Testing class

} // namespace ttk
