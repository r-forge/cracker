pb.modified <-
function (title = "R progress bar", label = "", min = 0, max = 1, 
    initial = 0, width = 300,.col = "#d9d9d9",.family = "helvetica") 
{
    useText <- FALSE
    have_ttk <- as.character(tcl("info", "tclversion")) >= "8.5"
    if (!have_ttk && as.character(tclRequire("PBar")) == "FALSE") 
        useText <- TRUE
    .win <- tktoplevel(background = .col)
    .val <- initial
    .killed <- FALSE
    tkwm.geometry(.win, sprintf("%dx75", width + 40))
    tkwm.resizable(.win, "FALSE","FALSE")

    tkwm.title(.win, title)
    fn <- tkfont.create(family = .family, size = 12)
    if (useText) {
        .lab <- tklabel(.win, text = label, font = fn, padx = 20)
        tkpack(.lab, side = "left")
        fn2 <- tkfont.create(family = ,family, size = 12)
        .vlab <- tklabel(.win, text = "0%", font = fn2, padx = 20)
        tkpack(.vlab, side = "right")
        up <- function(value) {
            if (!is.finite(value) || value < min || value > max) 
                return()
            .val <<- value
            tkconfigure(.vlab, text = sprintf("%d%%", round(100 * 
                (value - min)/(max - min))))
        }
    }
    else {
        .lab <- tklabel(.win, text = label, font = fn, pady = 10,background = .col)
        .tkval <- tclVar(0)
        tkpack(.lab, side = "top")
        tkpack(tklabel(.win, text = "", font = fn,background = .col), side = "bottom")
        pBar <- if (have_ttk) 
            ttkprogressbar(.win, length = width, variable = .tkval)
        else tkwidget(.win, "ProgressBar", width = width, variable = .tkval)
        tkpack(pBar, side = "bottom")
        up <- function(value) {
            if (!is.finite(value) || value < min || value > max) 
                return()
            .val <<- value
            tclvalue(.tkval) <<- 100 * (value - min)/(max - min)
        }
    }
    getVal <- function() .val
    kill <- function() if (!.killed) {
        tkdestroy(.win)
        .killed <<- TRUE
    }
    title <- function(title) tkwm.title(.win, title)
    lab <- function(label) tkconfigure(.lab, text = label)
    tkbind(.win, "<Destroy>", kill)
    up(initial)
    structure(list(getVal = getVal, up = up, title = title, label = lab, 
        kill = kill), class = "tkProgressBar")
}
