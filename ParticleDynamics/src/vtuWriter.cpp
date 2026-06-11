#include <fstream>
#include <string>

void writeVTU(const std::string& filename, const double* pos, int N, double radius)
{
    std::ofstream f(filename);
    f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
    f << "  <UnstructuredGrid>\n";
    f << "  <Piece NumberOfPoints=\"" << N << "\" NumberOfCells=\"0\">\n";
    f << "    <PointData Scalars=\"Radius\">\n";
    f << "      <DataArray type=\"Float64\" Name=\"Radius\" format=\"ascii\">\n";
    for (int i = 0; i < N; i++)
        f << "        " << radius << "\n";
    f << "      </DataArray>\n";
    f << "    </PointData>\n";
    f << "    <Points>\n";
    f << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < N; i++)
        f << "        " << pos[3*i+0] << " " << pos[3*i+1] << " " << pos[3*i+2] << "\n";
        // f << "        " << pos[3*i+0] << " " << pos[3*i+1] << " " << 0 << "\n";
    f << "      </DataArray>\n";
    f << "    </Points>\n";
    f << "    <Cells>\n";
    f << "      <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"/>\n";
    f << "      <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"/>\n";
    f << "      <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"/>\n";
    f << "    </Cells>\n";
    f << "  </Piece>\n";
    f << "  </UnstructuredGrid>\n";
    f << "</VTKFile>\n";
}

