#' Run the default model using UK population data from polymod and
#' demography groups (0, 20, 65)
test_data_default_vaccinated <- function() {
    # load contact and population data from socialmixr::polymod
    polymod <- socialmixr::polymod
    contact_data <- socialmixr::contact_matrix(
        polymod,
        countries = "United Kingdom",
        age.limits = c(0, 20, 65),
        symmetric = TRUE
    )

    # prepare contact matrix
    contact_matrix <- t(contact_data$matrix)

    # prepare the demography vector
    demography_vector <- contact_data$demography$population
    names(demography_vector) <- rownames(contact_matrix)

    # initial conditions
    initial_i <- 1e-6
    initial_conditions <- c(
        S = 1 - initial_i, E = 0, I = initial_i, R = 0, V = 0
    )

    # build for all age groups
    initial_conditions <- rbind(
        initial_conditions,
        initial_conditions,
        initial_conditions
    )

    # assign rownames for clarity
    rownames(initial_conditions) <- rownames(contact_matrix)

    uk_population <- population(
        name = "UK",
        contact_matrix = contact_matrix,
        demography_vector = demography_vector,
        initial_conditions = initial_conditions
    )

    # prepare a vaccination object
    vaccinate_elders <- vaccination(
        name = "vaccinate elders",
        time_begin = matrix(100, nrow(contact_matrix)),
        time_end = matrix(250, nrow(contact_matrix)),
        nu = matrix(c(0, 0, 0.0001))
    )

    # run an epidemic model using `epidemic`
    output <- model_default_cpp(
        population = uk_population,
        vaccination = vaccinate_elders,
        time_end = 600, increment = 1.0
    )
}

test_that("Get the number of vaccinated and vaccination rate rescale factor for the default model", {
    # run epidemic model, expect no condition
    data <- test_data_default_vaccinated()
    expect_no_condition(
        vaccinated(data)
    )

    out <- vaccinated(data)


    # check for output type and contents
    expect_s3_class(out, "data.frame")


    expect_length(out, 4L)
    expect_named(
        out, c("demography_group", "vaccinated", "total", "rescale"),
        ignore.order = TRUE
    )

    # check for all positive values within the range 0 and total population size
    expect_true(
        all(
            out[["vaccinated"]] >= 0 & out[["vaccinated"]] <= out[["total"]] & out[["rescale"]] >= 1
        )
    )
})
