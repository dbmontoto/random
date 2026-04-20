#ifndef UTILITIES_H
#define UTILITIES_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>

namespace pcsaft {

/**
 * @brief Table class for formatted output
 */
class Table {
public:
    Table() = default;

    // Set table properties
    void setTitle(const std::string& title) { title_ = title; }
    void setHeaders(const std::vector<std::string>& headers) { headers_ = headers; }
    void setPrecision(int precision) { precision_ = precision; }
    void setColumnWidth(int width) { col_width_ = width; }

    // Add data rows
    void addRow(const std::vector<double>& row);
    void addRow(const std::vector<std::string>& row);
    void addMixedRow(const std::vector<std::string>& row); // For mixed string/number rows

    // Clear table
    void clear();

    // Print to console
    void print() const;

    // Export to CSV
    bool exportCSV(const std::string& filename) const;

    // Get as string
    std::string toString() const;

private:
    std::string title_;
    std::vector<std::string> headers_;
    std::vector<std::vector<std::string>> rows_;
    int precision_ = 6;
    int col_width_ = 15;

    void printSeparator(int total_width) const;
};

/**
 * @brief Utility functions
 */
namespace utils {

    /**
     * @brief Convert temperature units
     */
    double K_to_C(double T);
    double C_to_K(double T);
    double K_to_F(double T);
    double F_to_K(double T);

    /**
     * @brief Convert pressure units
     */
    double Pa_to_bar(double P);
    double bar_to_Pa(double P);
    double Pa_to_psi(double P);
    double psi_to_Pa(double P);
    double Pa_to_atm(double P);
    double atm_to_Pa(double P);

    /**
     * @brief Convert density units
     */
    double molm3_to_kgm3(double rho, double MW);  // [mol/m³] to [kg/m³]
    double kgm3_to_molm3(double rho, double MW);   // [kg/m³] to [mol/m³]

    /**
     * @brief Format double with specified precision
     */
    std::string formatDouble(double value, int precision = 6);

    /**
     * @brief Format scientific notation
     */
    std::string formatScientific(double value, int precision = 3);

    /**
     * @brief Create timestamp string
     */
    std::string getTimestamp();

    /**
     * @brief File operations
     */
    bool fileExists(const std::string& filename);
    std::string getFileExtension(const std::string& filename);

} // namespace utils

} // namespace pcsaft

#endif // UTILITIES_H
