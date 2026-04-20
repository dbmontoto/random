#include "thermo/pcsaft/utilities/utilities.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <sys/stat.h>

namespace pcsaft {

// ============================================================================
// Table Implementation
// ============================================================================

void Table::addRow(const std::vector<double>& row) {
    std::vector<std::string> str_row;
    str_row.reserve(row.size());

    for (double val : row) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(precision_) << val;
        str_row.push_back(oss.str());
    }

    rows_.push_back(str_row);
}

void Table::addRow(const std::vector<std::string>& row) {
    rows_.push_back(row);
}

void Table::addMixedRow(const std::vector<std::string>& row) {
    rows_.push_back(row);
}

void Table::clear() {
    title_.clear();
    headers_.clear();
    rows_.clear();
}

void Table::printSeparator(int total_width) const {
    std::cout << std::string(total_width, '-') << "\n";
}

void Table::print() const {
    int num_cols = headers_.size();
    int total_width = col_width_ * num_cols + (num_cols - 1);

    // Print title
    if (!title_.empty()) {
        std::cout << "\n";
        printSeparator(total_width);
        std::cout << title_ << "\n";
        printSeparator(total_width);
    }

    // Print headers
    if (!headers_.empty()) {
        for (size_t i = 0; i < headers_.size(); ++i) {
            std::cout << std::setw(col_width_) << std::right << headers_[i];
            if (i < headers_.size() - 1) std::cout << " ";
        }
        std::cout << "\n";
        printSeparator(total_width);
    }

    // Print rows
    for (const auto& row : rows_) {
        for (size_t i = 0; i < row.size() && i < headers_.size(); ++i) {
            std::cout << std::setw(col_width_) << std::right << row[i];
            if (i < row.size() - 1) std::cout << " ";
        }
        std::cout << "\n";
    }

    if (!rows_.empty()) {
        printSeparator(total_width);
    }
    std::cout << "\n";
}

bool Table::exportCSV(const std::string& filename) const {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << " for writing\n";
        return false;
    }

    // Write title as comment
    if (!title_.empty()) {
        file << "# " << title_ << "\n";
    }

    // Write headers
    for (size_t i = 0; i < headers_.size(); ++i) {
        file << headers_[i];
        if (i < headers_.size() - 1) file << ",";
    }
    file << "\n";

    // Write rows
    for (const auto& row : rows_) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i < row.size() - 1) file << ",";
        }
        file << "\n";
    }

    file.close();
    return true;
}

std::string Table::toString() const {
    std::ostringstream oss;
    int num_cols = headers_.size();
    int total_width = col_width_ * num_cols + (num_cols - 1);

    // Title
    if (!title_.empty()) {
        oss << "\n" << std::string(total_width, '-') << "\n";
        oss << title_ << "\n";
        oss << std::string(total_width, '-') << "\n";
    }

    // Headers
    if (!headers_.empty()) {
        for (size_t i = 0; i < headers_.size(); ++i) {
            oss << std::setw(col_width_) << std::right << headers_[i];
            if (i < headers_.size() - 1) oss << " ";
        }
        oss << "\n" << std::string(total_width, '-') << "\n";
    }

    // Rows
    for (const auto& row : rows_) {
        for (size_t i = 0; i < row.size() && i < headers_.size(); ++i) {
            oss << std::setw(col_width_) << std::right << row[i];
            if (i < row.size() - 1) oss << " ";
        }
        oss << "\n";
    }

    if (!rows_.empty()) {
        oss << std::string(total_width, '-') << "\n\n";
    }

    return oss.str();
}

// ============================================================================
// Utility Functions
// ============================================================================

namespace utils {

// Temperature conversions
double K_to_C(double T) {
    return T - 273.15;
}

double C_to_K(double T) {
    return T + 273.15;
}

double K_to_F(double T) {
    return (T - 273.15) * 9.0 / 5.0 + 32.0;
}

double F_to_K(double T) {
    return (T - 32.0) * 5.0 / 9.0 + 273.15;
}

// Pressure conversions
double Pa_to_bar(double P) {
    return P / 1e5;
}

double bar_to_Pa(double P) {
    return P * 1e5;
}

double Pa_to_psi(double P) {
    return P / 6894.76;
}

double psi_to_Pa(double P) {
    return P * 6894.76;
}

double Pa_to_atm(double P) {
    return P / 101325.0;
}

double atm_to_Pa(double P) {
    return P * 101325.0;
}

// Density conversions
double molm3_to_kgm3(double rho, double MW) {
    return rho * MW / 1000.0; // MW in g/mol, convert to kg/mol
}

double kgm3_to_molm3(double rho, double MW) {
    return rho * 1000.0 / MW;
}

// Formatting
std::string formatDouble(double value, int precision) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << value;
    return oss.str();
}

std::string formatScientific(double value, int precision) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(precision) << value;
    return oss.str();
}

std::string getTimestamp() {
    std::time_t now = std::time(nullptr);
    char buf[100];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
    return std::string(buf);
}

// File operations
bool fileExists(const std::string& filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

std::string getFileExtension(const std::string& filename) {
    size_t pos = filename.find_last_of('.');
    if (pos == std::string::npos) {
        return "";
    }
    return filename.substr(pos + 1);
}

} // namespace utils

} // namespace pcsaft
