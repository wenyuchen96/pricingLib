#include <catch2/catch_all.hpp>
#include <catch2/catch_session.hpp>

int main(int argc, char *argv[])
{
    Catch::Session session;
    return session.run(argc, argv);
}
