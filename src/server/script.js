function change_dir() {
    var dir_path_input = document.getElementById("dir_path");
    document.getElementById("dir_listing").src = "dir_listing?dir_path=" + dir_path_input.value;
}

function generate_plots() {
    var inputs = document.getElementById("plot_command_pane").getElementsByTagName("input");

    var plot_command = "";
    for (var i = 0; i < inputs.length; i++) {
        plot_command += "--" + inputs[i].id + " " + inputs[i].value + " ";
    }
    console.log(plot_command);

    let xhr = new XMLHttpRequest();
    xhr.open("POST", "/generate_plots");
    xhr.setRequestHeader("Accept", "application/json");
    xhr.setRequestHeader("Content-Type", "application/json");

    xhr.onreadystatechange = function () {
        if (xhr.readyState === 4) {
            console.log(xhr.status);
            console.log(xhr.responseText);

            if (xhr.status == 200) {
                document.getElementById("dir_path").value = document.getElementById("output_dir").value;
                document.getElementById("list_files").click();

                document.getElementById("command_history").innerHTML += "<li>" + plot_command + "</li>"; 
            }
        }
    };

    xhr.send(plot_command);
}